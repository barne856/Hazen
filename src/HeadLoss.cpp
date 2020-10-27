#include "HeadLoss.hpp"
#include <algorithm>
#include <numeric>
#define _USE_MATH_DEFINES
#include <math.h>

namespace hazen {

double vertical_shaft_loss(HydraulicShape *shape, FrictionMethod *friction,
                           double up_inv, double dn_inv, double Q, double h,
                           std::vector<std::pair<double, double>> &HGL) {
  double A = shape->flow_area(shape->get_max_depth());
  double V = A <= 0.0 ? 0.0 : Q / A;
  double hv = V * V / (2.0 * g);
  double Sf = friction->friction_slope(shape, Q, shape->get_max_depth());

  // vertical pipe with discharge at the top, always calculated full length
  if (dn_inv > up_inv) {
    double P = shape->wetted_perimeter(shape->get_max_depth());
    Rectangle rect = Rectangle(P, 1.0, true);
    double d2 = std::max(critical_depth(&rect, Q), h - dn_inv);
    double L = dn_inv - up_inv;
    double H1 = dn_inv + d2 + hv + Sf * L;
    double d1 = H1 - up_inv - hv;
    // write HGL
    HGL.clear();
    HGL.push_back({0.0, dn_inv + d2});
    HGL.push_back({0.0, up_inv + d1});
    return H1;
  }
  // vertical pipe with discharge at the bottom
  else if (dn_inv < up_inv) {
    double depth_to_invert = h - dn_inv;
    if (depth_to_invert < 0.0) {
      // write HGL
      HGL.clear();
      HGL.push_back({0.0, h});
      HGL.push_back({0.0, h});
      return h + hv;
    }
    double h1 = h + Sf * (depth_to_invert) / (1.0 - Sf);
    // write HGL
    HGL.clear();
    HGL.push_back({0.0, h});
    HGL.push_back({0.0, h1});
    return h1 + hv;
  } else {
    // write HGL
    HGL.clear();
    HGL.push_back({0.0, h});
    HGL.push_back({0.0, h});
    return h + hv;
  }
}

double gvf_backwater_loss(HydraulicShape *shape, FrictionMethod *friction,
                          vec3 up_inv, vec3 dn_inv, double Q, double h,
                          double ds, double &jump_x,
                          std::vector<std::pair<double, double>> &HGL) {
  double S = slope(dn_inv, up_inv);     // bottom slope.
  double a = sqrt(S * S + 1);           // large angle correction
  double d2 = (h - dn_inv.z);           // downstream control point depth.
  double dc = critical_depth(shape, Q); // critical depth
  double dn = normal_depth(shape, friction, S, Q); // normal depth
  bool is_steep = dn < dc;
  if (d2 < dc && is_steep) {
    // cannot perform backwater calculations on supercritical flow
    jump_x = 0.0;
    return nan("");
  }
  // downstream control point depth cannot be less than critical depth.
  d2 = std::max(dc, d2);
  // // if normal depth is infinite, starting depth is at the crown.
  // if (dn == std::numeric_limits<double>::infinity()) {
  //   d2 = a * shape->get_max_depth();
  // }
  // change in depth with change in distance along passage.
  std::function<double(double, double)> dddx =
      [shape, friction, Q, a, S, ds, dn, dc, is_steep](double x,
                                                       double h) -> double {
    double d = h / a; // large angle correction
    double Sf = friction->friction_slope(shape, Q, d);
    double Fr = shape->froude(Q, d);
    if (abs(1.0 - Fr * Fr) < ds || d < dc) {
      // remove discontinuity around critical depth
      Fr = 0.0;
    }
    double test_slope = (Sf - S) / (1.0 - Fr * Fr);
    // if new depth will be below normal depth, jump to normal depth
    if (d + test_slope * ds < dn && is_steep) {
      return (dn - d) / ds;
    }
    return test_slope;
  };
  HGL.clear();
  HGL.push_back({0.0, dn_inv.z + d2});
  // integrate
  double d = d2;
  double s = 0.0;
  double L = length(dn_inv, up_inv);
  while (s + ds <= L) {
    double x = s / a;
    double dx = ds / a;
    // calculate next h using Runge-Kutta method.
    d = RK4(dddx, d, x, dx);
    // increment x by one step.
    s += ds;
    x = s / a;
    HGL.push_back({x, dn_inv.z + S * x + d * a});
    if (d <= dc && is_steep) {
      // hydraulic jump has occured.
      jump_x = x;
      return nan("");
    }
  }
  // check if there is a step less than dx left at the end.
  if (s < L) {
    // compute a partial step to reach the end of the passage.
    ds = L - s;
    double x = s / a;
    double dx = ds / a;
    // calculate next h using Runge-Kutta method.
    d = RK4(dddx, d, x, dx);
    // increment x by one step.
    s += ds;
    x = s / a;
    HGL.push_back({x, dn_inv.z + S * x + d * a});
  }
  jump_x = nan("");
  // return the upstream energy head.
  double V1 = Q / shape->flow_area(d);
  double hv = V1 * V1 / (2.0 * g);
  return d + up_inv.z + hv;
}

double gvf_frontwater_loss(HydraulicShape *shape, FrictionMethod *friction,
                           vec3 up_inv, vec3 dn_inv, double Q, double h,
                           double ds, double jump_x,
                           std::vector<std::pair<double, double>> &HGL) {
  double dc = critical_depth(shape, Q);
  double d1 = std::max(0.0, std::min(dc, h - up_inv.z));
  double A = shape->flow_area(d1);
  double V1 = A <= 0.0 ? 0.0 : Q / A;
  double H1 = d1 + V1 * V1 / (2.0 * g) + up_inv.z;
  jump_x = isnan(jump_x) ? 0.0 : jump_x;
  double S = slope(dn_inv, up_inv);                // bottom slope.
  double a = sqrt(S * S + 1);                      // large angle correction
  double dn = normal_depth(shape, friction, S, Q); // normal depth

  std::function<double(double, double)> dddx =
      [shape, friction, Q, a, S, ds, dc](double x, double h) -> double {
    // double d = h / a; // large angle correction
    // double Sf = friction->friction_slope(shape, Q, d);
    // // approximate without froude number since there is no hydraulic jump
    // // possible.
    // return (Sf - S);

    double d = h / a; // large angle correction
    double Sf = friction->friction_slope(shape, Q, d);
    double Fr = shape->froude(Q, d);
    if (abs(1.0 - Fr * Fr) < ds || d > dc) {
      // remove discontinuity around critical depth
      Fr = 2.0;
    }
    return (S - Sf) / (1.0 - Fr * Fr);
  };
  // integrate
  double d = d1;
  double s = 0.0;
  double L = length(dn_inv, up_inv);
  HGL.clear();
  HGL.push_back({L / a, up_inv.z + d});
  L = L - jump_x * a;
  while (s + ds <= L) {
    double x = (L - s) / a + jump_x;
    double dx = ds / a;
    // calculate next h using Runge-Kutta method.
    d = RK4(dddx, d, x, dx);
    // increment x by one step.
    s += ds;
    x = (L - s) / a + jump_x;
    HGL.push_back({x, dn_inv.z + S * x + d * a});
  }
  // check if there is a step less than dx left at the end.
  if (s < L) {
    // compute a partial step to reach the end of the passage.
    ds = L - s;
    double x = (L - s) / a + jump_x;
    double dx = ds / a;
    // calculate next h using Runge-Kutta method.
    d = RK4(dddx, d, x, dx);
    // increment x by one step.
    s += ds;
    x = (L - s) / a + jump_x;
    HGL.push_back({x, dn_inv.z + S * x + d * a});
  }
  // return the upstream energy head.
  return H1;
}

double opening_loss(HydraulicShape *shape, double Cd, double opening_invert,
                    double Q, double h, double dy,
                    std::vector<std::pair<double, double>> &HGL,
                    double percent_open) {
  // get constants
  double D = shape->get_max_depth() * percent_open;
  double d2 = h - opening_invert;
  // define velocity as function of y (distance up from opening invert)
  std::function<double(double, double)> V = [=](double y, double d1) -> double {
    if (0.0 <= y && y < std::max(0.0, d2)) {
      return sqrt(2.0 * g * (d1 - d2));
    } else if (std::max(0.0, d2) <= y && y <= std::min(d1, D)) {
      return sqrt(2.0 * g * (d1 - y));
    } else {
      return 0.0;
    }
  };
  // define flow change as function of y (distance up from opening invert)
  std::function<double(double, double, double)> dQdy =
      [=](double y, double Q, double d1) -> double {
    return Cd * V(y, d1) * shape->top_width(y);
  };
  // define flow as a function of d1
  std::function<double(double)> Q_d1 = [=](double d1) -> double {
    return integrate(
        RK4, std::bind(dQdy, std::placeholders::_1, std::placeholders::_2, d1),
        0.0, 0.0, std::min(d1, D), dy);
  };
  // goal seek Q_d1 for d1
  double dc = critical_depth(shape, Q);
  double a = std::max(dc, d2);
  double b = a + 0.1;
  double d1 = find_goal_secant(Q, a, b, 0.00000001, Q_d1);
  // write HGL
  HGL.clear();
  HGL.push_back({0.0, opening_invert + d2});
  HGL.push_back({0.0, opening_invert + d1});
  // return total energy head upstream of opening "H1"
  return d1 + opening_invert;
}

double transition_loss(HydraulicShape *up_shape, HydraulicShape *dn_shape,
                       double up_inv, double dn_inv, double Q, double K,
                       double h, std::vector<std::pair<double, double>> &HGL) {
  // section 1 is upstream, section 2 is downstream.
  double A2 = dn_shape->flow_area(h - dn_inv);
  double V2 = A2 <= 0.0 ? 0.0 : Q / A2;
  double H2 = h + V2 * V2 / (2.0 * g);
  // Check that downstream depth is subcritical. If not, set downstream depth
  // to the critical depth.
  double d2 = std::max(h - dn_inv, critical_depth(dn_shape, Q));
  V2 = Q / dn_shape->flow_area(d2);

  // objective function to solve for d1.
  std::function<double(double)> objective = [=](double d1) {
    double V1 = Q / up_shape->flow_area(d1);
    if (V1 < V2) {
      // use velocity in the vena contracta
      double mu = 0.63 + 0.37 * pow(V1 / V2, 3.0); // according to Weisbach
      double V3 = V2 / mu;
      return dn_inv + d2 + V2 * V2 / (2.0 * g) +
             K * (V3 - V2) * (V3 - V2) / (2.0 * g) - up_inv - d1 -
             V1 * V1 / (2.0 * g);
    }
    return dn_inv + d2 + V2 * V2 / (2.0 * g) +
           K * (V1 - V2) * (V1 - V2) / (2.0 * g) - up_inv - d1 -
           V1 * V1 / (2.0 * g);
  };

  double d1 = find_goal_secant(0.0, d2, d2 + 0.1, 0.0000001, objective);
  // if no solution, set d1 to the critical depth.
  double d1_critical = critical_depth(up_shape, Q);
  d1 = isnan(d1) || d1 < d1_critical ? d1_critical : d1;
  // write results to HGL vector
  HGL.clear();
  HGL.push_back({0.0, dn_inv + d2});
  HGL.push_back({0.0, up_inv + d1});
  // return total energy head upstream of transition
  double V1 = Q / up_shape->flow_area(d1);
  double H1 = up_inv + d1 + V1 * V1 / (2.0 * g);
  return H1;
}

double benching_coefficient(BENCH_CONFIGURATION bench, double D, double d) {
  double x0 = 1.0, x1 = 2.5, y0, y1;
  switch (bench) {
  case BENCH_CONFIGURATION::FLAT:
    y0 = -0.05;
    y1 = -0.05;
    break;
  case BENCH_CONFIGURATION::DEPRESSED:
    y0 = 0.00;
    y1 = 0.00;
    break;
  case BENCH_CONFIGURATION::HALF:
    y0 = -0.85;
    y1 = -0.05;
    break;
  case BENCH_CONFIGURATION::FULL:
    y0 = -0.93;
    y1 = -0.25;
    break;
  case BENCH_CONFIGURATION::IMPROVED:
    y0 = -0.98;
    y1 = -0.60;
    break;
  default:
    y0 = 0.00;
    y1 = 0.00;
    break;
  }
  if (d <= 0) {
    d = D;
  }
  double x;
  if (D <= 0.0) {
    x = 0.0;
  } else {
    x = D / d;
  }

  // If D/d >= 2.5, submerged
  if (x >= 2.5) {
    return y1;
  }
  // If D/d <= 1.0, unsubmerged
  if (x <= 1.0) {
    return y0;
  }
  // linearly interpolate y
  double a = (x - x0) / (x1 - x0);
  return y0 * (1.0 - a) + y1 * a;
}

double manhole_loss(double E_ai, double E_i, double Q_0, double D_0,
                    double invert_0, BENCH_CONFIGURATION bench,
                    std::vector<double> Q, std::vector<double> invert,
                    std::vector<double> theta,
                    std::vector<std::pair<double, double>> &HGL) {
  // Error handling
  if (Q_0 <= 0.0 || D_0 <= 0.0) {
    return E_ai;
  }

  // find all plunging flows j and non plunging flows k
  std::vector<double> theta_j{};
  std::vector<double> Q_j{};
  std::vector<double> invert_k{};
  std::vector<double> Q_k{};
  for (int i = 0; i < invert.size(); i++) {
    if (invert[i] < E_ai) {
      // pipe is submerged
      theta_j.push_back(theta[i]);
      Q_j.push_back(Q[i]);
    } else {
      // pipe is not submerged
      invert_k.push_back(invert[i]);
      Q_k.push_back(Q[i]);
    }
  }

  // Benching ------------------------------------------------------------------
  double C_B = benching_coefficient(bench, D_0, E_ai - invert_0);

  // Angled Inflow -------------------------------------------------------------
  double theta_w = 180.0;
  double Q_j_sum = std::accumulate(std::begin(Q_j), std::end(Q_j), 0.0);
  if (Q_j_sum > 0.0) {
    theta_w = std::inner_product(std::begin(Q_j), std::end(Q_j),
                                 std::begin(theta_j), 0.0) /
              Q_j_sum;
  }
  double C_theta = 4.5 * (Q_j_sum / Q_0) * cos((M_PI / 180.0) * theta_w / 2.0);

  // Plunging ------------------------------------------------------------------
  std::vector<double> h_k{};
  for (int k = 0; k < invert_k.size(); k++) {
    double h = (invert_k[k] - E_ai) / D_0;
    if (h > 10.0 * D_0) {
      h = 10.0 * D_0;
    }
    h_k.push_back(h);
  }
  double C_P =
      std::inner_product(std::begin(Q_k), std::end(Q_k), std::begin(h_k), 0.0) /
      Q_0;

  // Return results
  double H_a = std::max(0.0, (C_B + C_theta + C_P) * (E_ai - E_i));
  HGL.clear();
  HGL.push_back({0.0, E_ai});
  HGL.push_back({0.0, E_ai + H_a});
  return E_ai + H_a;
}


double bar_screen_loss(HydraulicShape *upstream_channel_shape,
                       double channel_invert,
                       HydraulicShape *bar_screen_opening_shape, int n,
                       double Q, double h) {
  std::function<double(double)> bar_screen_approach_velocity = [=](double hl) {
    double upstream_wse = h + hl;
    double upstream_depth = upstream_wse - channel_invert;
    double A = upstream_channel_shape->flow_area(upstream_depth);
    return Q / A;
  };
  std::function<double(double)> objective = [=](double hl) {
    double u = bar_screen_approach_velocity(hl);
    double d = h-channel_invert; // depth of flow though the bar screen.
    double q = Q/n;
    double A = bar_screen_opening_shape->flow_area(d);
    double v = q/A;
    return (v * v - u * u) / (0.7 * 2 * g) - hl;
  };
  double hl = find_goal_secant(0.0, 1.0, 0.5, 0.0000001, objective);
  double u = bar_screen_approach_velocity(hl);
  double H = h+hl+(u*u)/(2.0*g); // upstream energy head in the channel.
  return H;
}

} // namespace hazen