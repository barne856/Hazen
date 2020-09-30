#include "HeadLoss.hpp"
#include <algorithm>

namespace hazen {

double gvf_backwater_loss(HydraulicShape *shape, FrictionMethod *friction,
                          vec3 up_inv, vec3 dn_inv, double Q, double H,
                          double dx, double &jump_x, double x_start,
                          double x_final) {
  double x = x_start;               // starting distance in the x direction
  double S = slope(dn_inv, up_inv); // bottom slope.
  if (Q <= 0.0) {
    jump_x = nan("");
    return std::max(S * x + dn_inv.z, H);
  }
  if (S == std::numeric_limits<double>::infinity()) {
    // vertical pipe with discharge at the bottom, always calculates full length
    double Sf = friction->friction_slope(shape, Q, shape->get_max_depth());
    double depth_to_invert = H - dn_inv.z;
    if (depth_to_invert < 0.0) {
      return H;
    }
    return H + Sf * (depth_to_invert) / (1.0 - Sf);
  }
  if (S == -std::numeric_limits<double>::infinity()) {
    // vertical pipe with discharge at the top, always calculated full length
    double P = shape->wetted_perimeter(shape->get_max_depth());
    Rectangle rect = Rectangle(P, 100.0);
    double dc = critical_depth(&rect, Q);
    H = std::max(dn_inv.z + dc, H);
    double Sf = friction->friction_slope(shape, Q, shape->get_max_depth());
    double L = length(dn_inv, up_inv);
    return H + Sf * L;
  }
  double a = sqrt(S * S + 1);           // large angle correction
  double h = (H - (S * x + dn_inv.z));  // downstream control point depth.
  double dc = critical_depth(shape, Q); // critical depth
  double dn = normal_depth(shape, friction, S, Q); // normal depth
  bool is_steep = dn < dc;
  if (h < dc && is_steep) {
    // cannot perform backwater calculations on supercritical flow
    jump_x = x;
    return H;
  }
  // downstream control point depth cannot be less than critical depth.
  h = std::max(dc, h);
  if (dn == std::numeric_limits<double>::infinity()) {
    h = a * shape->get_max_depth();
  }
  double L = horizontal_length(dn_inv, up_inv);
  L = std::min(x_final, L);
  std::function<double(double, double)> dhdx =
      [shape, friction, Q, a, S, dx, dn, dc, is_steep](double x,
                                                       double h) -> double {
    double d = h / a; // large angle correction
    double Sf = friction->friction_slope(shape, Q, d);
    double Fr = shape->froude(Q, d);
    if (abs(1.0 - Fr * Fr) < 0.01 * dx / d || d < dc) {
      // remove discontinuity around critical depth
      Fr = 0.0;
    }
    double test_slope = (Sf - S) / (1.0 - Fr * Fr);
    // if new depth will be below normal depth, jump to normal depth
    if (test_slope * dx / a > d - dn && is_steep) {
      return a * (d - dn) / dx;
    }
    return test_slope;
  };
  double d;
  while (x <= L - dx) {
    // calculate next h using Runge-Kutta method.
    h = RK4(dhdx, h, x, dx);
    // increment x by one step.
    x += dx;
    d = h / a; // flow depth normal to slope.
    if (d < dc && is_steep) {
      // hydraulic jump has occured.
      jump_x = x;
      // return HGL of normal depth at jump location
      return a * dn + S * x + dn_inv.z;
    }
  }
  // check if there is a step less than dx left at the end.
  if (x < L) {
    // compute a partial step to reach the end of the passage.
    dx = L - x;
    h = RK4(dhdx, h, x, dx);
    x += dx;
  }
  jump_x = nan("");
  // return the upstream water surface elevation.
  return h + S * x + dn_inv.z;
}

double gvf_frontwater_loss(HydraulicShape *shape, FrictionMethod *friction,
                           vec3 up_inv, vec3 dn_inv, double Q, double H,
                           double dx, double x_start, double x_final) {
  double x = x_start;               // starting distance in the x direction
  double S = slope(dn_inv, up_inv); // bottom slope.
  if (Q <= 0.0) {
    return std::max(up_inv.z - S * x, H);
  }
  if (S == std::numeric_limits<double>::infinity()) {
    // vertical pipe with discharge at the bottom, always calculates full length
    double Sf = friction->friction_slope(shape, Q, shape->get_max_depth());
    double depth_to_invert = H - dn_inv.z;
    if (depth_to_invert < 0.0) {
      return H;
    }
    return H + Sf * (depth_to_invert) / (1.0 - Sf);
  }
  if (S == -std::numeric_limits<double>::infinity()) {
    // vertical pipe with discharge at the top, always calculated full length
    double P = shape->wetted_perimeter(shape->get_max_depth());
    Rectangle rect = Rectangle(P, 100.0);
    double dc = critical_depth(&rect, Q);
    H = std::max(dn_inv.z + dc, H);
    double Sf = friction->friction_slope(shape, Q, shape->get_max_depth());
    double L = length(dn_inv, up_inv);
    return H + Sf * L;
  }
  double a = sqrt(S * S + 1);                      // large angle correction
  double dc = critical_depth(shape, Q);            // critical depth
  double dn = normal_depth(shape, friction, S, Q); // normal depth
  // upstream control point depth cannot be more than critical depth or less
  // than 0.0.
  double h = (H - (up_inv.z - S * x));
  h = std::min(h, dc);
  if (h <= 0.0) {
    h = dc;
  }
  double L = horizontal_length(dn_inv, up_inv);
  L = std::min(x_final, L);
  std::function<double(double, double)> dhdx =
      [shape, friction, Q, a, S, dx, dc](double x, double h) -> double {
    double d = h / a; // large angle correction
    double Sf = friction->friction_slope(shape, Q, d);
    // approximate without froude number
    return (Sf - S);

    // double d = h / a; // large angle correction
    // double Sf = friction->friction_slope(shape, Q, d);
    // double Fr = shape->froude(Q, d);
    // if (abs(1.0 - Fr * Fr) < 0.01 * dx / d || d > dc) {
    //   // remove discontinuity around critical depth
    //   Fr = 2.0;
    // }
    // return (S - Sf) / (1.0 - Fr * Fr);
  };
  while (x <= L - dx) {
    // calculate next h using Runge-Kutta method.
    h = RK4(dhdx, h, x, dx);
    // increment x by one step.
    x += dx;
  }
  // check if there is a step less than dx left at the end.
  if (x < L) {
    // compute a partial step to reach the end of the passage.
    dx = L - x;
    h = RK4(dhdx, h, x, dx);
    x += dx;
  }
  // return the downstream water surface elevation.
  return h + up_inv.z - S * x;
}
// approximates weir and orifice flow according to downstream conditions.
// works with free and submerged orifices or weirs
double opening_loss(HydraulicShape *shape, double Cd, double opening_invert,
                    double Q, double H, double dy, double percent_open) {
  // downstream head over the opening invert.
  double h2 = std::max(H - opening_invert, 0.0);
  // function to compute Q_guess from upstream head
  std::function<double(double)> objective = [shape, Cd, h2, dy, percent_open,
                                             Q](double h1) -> double {
    if (h1 < h2) {
      return nan("");
    }
    // compute the amount of flow contributed from the orifice-type flow if h2 >
    // 0
    double q = sqrt(2 * g * (h1 - h2)) * shape->flow_area(h2);
    double h = h2;
    // function to approximate the derivative of flow
    std::function<double(double, double)> dQdh =
        [shape, h1](double h, double x) -> double {
      if (h1 < h) {
        return nan("");
      }
      return sqrt(2 * g * (h1 - h)) * shape->top_width(h);
    };
    double h_end = std::min(h1, shape->get_max_depth() * percent_open);
    while (h + dy <= h_end) {
      q = RK4(dQdh, q, h, dy);
      h += dy;
    }
    if (h < h_end) {
      double dh = h_end - h;
      q = RK4(dQdh, q, h, dh);
      h += dh;
    }
    return q;
  };
  double b = 2.0 * (h2 + 1.0);
  while (objective(b) < Q) {
    b *= 2.0;
  }
  double h1 = find_goal_bisection(Q, h2, b, 0.0000001, objective);
  return opening_invert + h1;
}

double manhole_loss(double E_ai, double E_i, double manhole_invert,
                    double Q_total, double D_outlet, double C_B,
                    std::vector<double> Q_j, std::vector<double> invert_j,
                    std::vector<double> theta_j) {
  double H_B = C_B * (E_ai - E_i); // headloss due to benching
  // find all non-plunging flows
  std::vector<double> non_plunging_invert_j{};
  std::vector<double> non_plunging_theta_j{};
  std::vector<double> non_plunging_Q_j{};
  for (int i = 0; i < invert_j.size(); i++) {
    if (invert_j[i] < E_ai) {
      // pipe is submerged
      non_plunging_invert_j.push_back(invert_j[i]);
      non_plunging_theta_j.push_back(theta_j[i]);
      non_plunging_Q_j.push_back(Q_j[i]);
    }
  }
  // theta_w
  double numerator;
  for (int i = 0; i < non_plunging_theta_j.size(); i++) {
    numerator += non_plunging_Q_j[i] * non_plunging_theta_j[i];
  }
  double sum_of_non_plunging_flows;
  for (int i = 0; i < non_plunging_theta_j.size(); i++) {
    sum_of_non_plunging_flows += non_plunging_Q_j[i];
  }
  double theta_w = 180.0;
  if (non_plunging_theta_j.size() > 0) {
    theta_w = numerator / sum_of_non_plunging_flows;
  }
  double C_theta =
      4.5 * (sum_of_non_plunging_flows / Q_total) * cos(theta_w / 2.0);
  double H_theta = C_theta * (E_ai - E_i); // headloss due to angle of inflow.

  // find all plunging flows
  std::vector<double> plunging_invert_j{};
  std::vector<double> plunging_theta_j{};
  std::vector<double> plunging_Q_j{};
  std::vector<double> h_k{};
  for (int i = 0; i < invert_j.size(); i++) {
    if (invert_j[i] >= E_ai) {
      // pipe is plunging
      plunging_invert_j.push_back(invert_j[i]);
      plunging_theta_j.push_back(theta_j[i]);
      plunging_Q_j.push_back(Q_j[i]);
      h_k.push_back((invert_j[i] - E_ai) / D_outlet);
    }
  }
  double C_P = 0.0;
  for (int i = 0; i < h_k.size(); i++) {
    C_P += plunging_Q_j[i] * h_k[i];
  }
  C_P /= Q_total;
  double H_P = C_P * (E_ai - E_i);  // headloss due to plunging.
  double H_a = H_B + H_theta + H_P; // total headloss due to the manhole.
  return H_a + E_ai; // return the total energy head in the manhole.
}

// double pump_loss(std::vector<std::pair<double, double>> flow_head_curve,
//                  double discharge_elevation, double Q, double E) {
//   double h = E - discharge_elevation;
//   h = -std::min(h, 0.0);
//   return interp_1D(flow_head_curve, Q)+h;
// }

} // namespace hazen