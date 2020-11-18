#include "HeadLoss.hpp"
#include "HydraulicLinks.hpp"
#include "HydraulicUtil.hpp"
#include <algorithm>
#include <numeric>

#define _USE_MATH_DEFINES
#include <math.h>

namespace hazen {

double vertical_drop_loss(PassageLink *link, HydraulicNode *node) {
  // clear recorded HGL profile and initialize data
  link->HGL.clear();
  double D = link->cross_section_shape->get_max_depth();
  double h = link->get_water_surface_subcritical(node);
  double inv = link->invert(node);
  double Sf = link->friction_slope(D);
  double d = h - inv;
  double dh = Sf * d / (1.0 - Sf);
  double l = link->length();
  double hv = link->velocity_head(D);

  if (inv > h) { // check if discharge is open to atmosphere
    dh = 0.0;
  } else if (dh + d > l) { // check if upstream head will be above the pipe
    dh = Sf * l;
  }

  // record results
  double H = h + dh + hv;
  link->antinode(node)->H = H;
  link->HGL.push_back({0.0, h});
  link->HGL.push_back({0.0, h + dh});
  return H;
}

double vertical_rise_loss(PassageLink *link, HydraulicNode *node) {
  // clear recorded HGL profile and initialize data
  link->HGL.clear();
  double D = link->cross_section_shape->get_max_depth();
  double h = link->get_water_surface_subcritical(node);
  double inv = link->invert(node);
  double Sf = link->friction_slope(D);
  double l = link->length();
  double hv = link->velocity_head(D);
  double P = link->cross_section_shape->wetted_perimeter(D);
  Rectangle rect = Rectangle(P, 1.0, true);
  double d = std::max(critical_depth(&rect, abs(link->Q)), h - inv);
  double H = inv + d + hv + Sf * l;

  // record results
  link->antinode(node)->H = H;
  link->HGL.clear();
  link->HGL.push_back({0.0, h});
  link->HGL.push_back({0.0, H - hv});
  return H;
}

bool flowing_into_supercritical(HydraulicNode *node) {
  for (auto &index : node->out_links()) {
    if (auto passage = dynamic_cast<PassageLink *>(node->links[index].get())) {
      if (passage->is_supercritical) {
        return true;
      }
    }
  }
  return false;
}

double backwater_loss(PassageLink *link, HydraulicNode *node) {
  // clear recorded HGL profile and initialize data
  link->HGL.clear();
  link->is_supercritical = false;
  double h;
  if (link->is_steep(node) && flowing_into_supercritical(node)) {
    h = link->get_water_surface_supercritical(node);
  } else {
    h = link->get_water_surface_subcritical(node);
  }
  double S = link->slope(node);
  double a = sqrt(S * S + 1.0);
  double L = link->horizontal_length();
  double d = (h - link->invert(node)) / a;
  double dc = link->critical_depth();
  // if starting depth is less than the critical depth and the slope is steep,
  // the flow is supercritical and the system is determined by an upstream
  // control point
  if (link->is_steep(node) && (isnan(d) || d <= dc)) {
    link->HGL.push_back({0.0, nan("")});
    return nan("");
  }
  // flow must be subcritical.
  if (isnan(d)) {
    d = dc;
  }
  d = std::max(dc, d);
  h = d * a + link->invert(node);
  std::function<bool(double, double)> subcritical = [=](double x, double h) {
    double z = S * x + link->invert(node);
    double d = (h - z) / a;
    if (d > dc) {
      return true;
    } else {
      return false;
    }
  };
  // integrate until the end is reached or the flow becomes supercritical and a
  // hydraulic jump occurs.
  link->HGL =
      RKF45(std::bind(&PassageLink::hydrualic_slope_subcritical, link,
                      std::placeholders::_1, std::placeholders::_2, node),
            h, {0.0, L}, subcritical);
  // get the last entry in the results and return the total energy head at the
  // upstream node.
  h = link->HGL[link->HGL.size() - 1].second;
  double x = link->HGL[link->HGL.size() - 1].first;
  double z = link->invert(link->antinode(node)) - S * (L - x);
  double hv = link->velocity_head(h - z);
  return h + hv;
}

double frontwater_loss(PassageLink *link, HydraulicNode *node, double jump_x) {
  // clear recorded HGL profile and initialize data
  link->HGL.clear();
  link->is_supercritical = true;
  double h = link->get_upstream_water_surface(node);
  double S = link->slope(node);
  double L = link->horizontal_length();
  double a = sqrt(S * S + 1.0);
  double d = (h - link->invert(link->antinode(node))) / a;
  double dc = link->critical_depth();
  // passage must be steep with starting depth less than or equal to critical
  // depth.
  if (!link->is_steep(node)) {
    link->HGL.push_back({L, nan("")});
    return nan("");
  }
  // flow must be supercritical.
  d = std::min(dc, d);
  h = d * a + link->invert(link->antinode(node));
  std::function<bool(double, double)> supercritical = [=](double x, double h) {
    double z = S * x + link->invert(node);
    double d = (h - z) / a;
    if (d > dc) {
      return false;
    } else {
      return true;
    }
  };
  // integrate until the end is reached or the flow becomes supercritical and a
  // hydraulic jump occurs.
  link->HGL =
      RKF45(std::bind(&PassageLink::hydrualic_slope_supercritical, link,
                      std::placeholders::_1, std::placeholders::_2, node),
            h, {L, jump_x}, supercritical);
  // get the first entry in the results and return the total energy head at the
  // upstream node.
  h = link->HGL[0].second;
  double z = link->invert(link->antinode(node));
  double hv = link->velocity_head(h - z);
  return h + hv;
}

double opening_loss(OpeningLink *link, HydraulicNode *node) {
  // clear recorded HGL profile and initialize data
  link->HGL.clear();
  double h = link->get_water_surface_subcritical(node);
  double d2 = std::max(h - link->elevation, 0.0);
  link->HGL.push_back({0.0, h});
  // goal seek Q for d1
  double dc = critical_depth(link->cross_section_shape.get(), link->Q);
  double a = std::max(dc, d2);
  double b = a + 0.1;
  double d1 = find_goal_secant(
      abs(link->Q), a, b, 0.00000001,
      std::bind(&OpeningLink::flow, link, std::placeholders::_1, d2));
  h = d1 + link->elevation;
  link->HGL.push_back({0.0, h});
  return h;
}

double transition_loss(TransitionLink *link, HydraulicNode *node) {
  // clear recorded HGL profile and initialize data
  link->HGL.clear();
  double h = link->get_water_surface_subcritical(node);
  double d2 = h - link->invert(node);
  HydraulicNode *up_node = link->antinode(node);
  double V2 = link->velocity(node, d2);
  double z2 = link->invert(node);
  double z1 = link->invert(up_node);

  // objective function to solve for d1.
  std::function<double(double)> objective = [=](double d1) {
    double V1 = link->velocity(up_node, d1);
    if (V1 < V2) {
      // use velocity in the vena contracta
      double mu = 0.63 + 0.37 * pow(V1 / V2, 3.0); // according to Weisbach
      double V3 = V2 / mu;
      return z2 + d2 + V2 * V2 / (2.0 * g) +
             link->K * (V3 - V2) * (V3 - V2) / (2.0 * g) - z1 - d1 -
             V1 * V1 / (2.0 * g);
    }
    return z2 + d2 + V2 * V2 / (2.0 * g) +
           link->K * (V1 - V2) * (V1 - V2) / (2.0 * g) - z1 - d1 -
           V1 * V1 / (2.0 * g);
  };

  double d1 = find_goal_secant(0.0, d2, d2 + 0.1, 0.0000001, objective);
  // if no solution, set d1 to critical depth
  d1 =
      isnan(d1) || d1 < d2 ? critical_depth(link->shape(up_node), link->Q) : d1;
  // write results to HGL vector
  link->HGL.push_back({0.0, h});
  link->HGL.push_back({0.0, z1 + d1});
  // return total energy head upstream of transition
  return z1 + d1 + link->velocity_head(up_node, d1);
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
    double d = h - channel_invert; // depth of flow through the bar screen.
    double q = Q / n;
    double A = bar_screen_opening_shape->flow_area(d);
    double v = q / A;
    return (v * v - u * u) / (0.7 * 2.0 * g) - hl;
  };
  double hl = find_goal_secant(0.0, 1.0, 0.5, 0.0000001, objective);
  double u = bar_screen_approach_velocity(hl);
  double H =
      h + hl + (u * u) / (2.0 * g); // upstream energy head in the channel.
  return H;
}

} // namespace hazen