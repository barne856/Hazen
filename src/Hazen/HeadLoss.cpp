#include "Hazen/HydraulicLinks.hpp"
#include <algorithm>
#include <iostream>
#include <numeric>

#include <math.h>

namespace hazen {

Length vertical_drop_loss(PassageLink *link, HydraulicNode *node) {
  // clear recorded HGL profile and initialize data
  link->HGL.clear();
  if (link->Q.val == 0.0) {
    auto up_node = link->antinode(node);
    Length H = std::max(node->H, link->invert(up_node));
    link->antinode(node)->H = H;
    link->HGL.push_back({Length(0.0), node->H});
    link->HGL.push_back({Length(0.0), up_node->H});
    return H;
  }
  Length D = link->cross_section_shape->get_max_depth();
  Length h = link->get_water_surface_subcritical(node);
  Length inv = link->invert(node);
  Angle Sf = link->friction_slope(D);
  Length d = h - inv;
  Length dh = Sf * d / (Dimensionless(1.0) - Sf);
  Length l = link->length();
  Length hv = link->velocity_head(D);

  if (inv > h) { // check if discharge is open to atmosphere
    dh = Length(0.0);
  } else if (dh + d > l) { // check if upstream head will be above the pipe
    dh = Sf * l;
  }

  // record results
  Length H = h + dh + hv;
  link->antinode(node)->H = H;
  link->HGL.push_back({Length(0.0), h});
  link->HGL.push_back({Length(0.0), h + dh});
  return H;
}

Length vertical_rise_loss(PassageLink *link, HydraulicNode *node) {
  // clear recorded HGL profile and initialize data
  link->HGL.clear();
  if (link->Q.val == 0.0) {
    auto up_node = link->antinode(node);
    Length H = std::max(node->H, link->invert(node));
    link->antinode(node)->H = H;
    link->HGL.push_back({Length(0.0), node->H});
    link->HGL.push_back({Length(0.0), up_node->H});
    return H;
  }
  Length D = link->cross_section_shape->get_max_depth();
  Length h = link->get_water_surface_subcritical(node);
  Length inv = link->invert(node);
  Angle Sf = link->friction_slope(D);
  Length l = link->length();
  Length hv = link->velocity_head(D);
  Length P = link->cross_section_shape->wetted_perimeter(D);
  Rectangle rect = Rectangle(P, Length(1.0), true);
  Length d = Length(
      std::max(critical_depth(&rect, abs(link->Q)).val, (h.val - inv.val)));
  Length H = inv + d + hv + Sf * l;

  // record results
  link->antinode(node)->H = H;
  link->HGL.clear();
  link->HGL.push_back({Length(0.0), h});
  link->HGL.push_back({Length(0.0), H - hv});
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

Length backwater_loss(PassageLink *link, HydraulicNode *node) {
  // clear recorded HGL profile and initialize data
  link->HGL.clear();
  link->is_supercritical = false;
  if (link->Q.val == 0.0) {
    auto up_node = link->antinode(node);
    Length H =
        std::max(std::max(node->H, link->invert(node)), link->invert(up_node));
    link->antinode(node)->H = H;
    link->HGL.push_back({Length(0.0), node->H});
    Length x = (node->H - link->invert(node)) / link->slope(node);
    if (x.val > 0.0 && x < link->horizontal_length()) {
      link->HGL.push_back({x, node->H});
    }
    link->HGL.push_back({Length(link->horizontal_length()), up_node->H});
    return H;
  }
  Length h;
  if (link->is_steep(node) && flowing_into_supercritical(node)) {
    h = link->get_water_surface_supercritical(node);
  } else {
    h = link->get_water_surface_subcritical(node);
  }
  Angle S = link->slope(node);
  Dimensionless a = sqrt(S * S + Dimensionless(1.0));
  Length L = link->horizontal_length();
  Length d = (h - link->invert(node)) / a;
  Length dc = link->critical_depth();
  // if starting depth is less than the critical depth and the slope is steep,
  // the flow is supercritical and the system is determined by an upstream
  // control point
  if (link->is_steep(node) && (isnan(d.val) || d <= dc)) {
    link->HGL.push_back({Length(0.0), Length(nan(""))});
    return Length(nan(""));
  }
  // flow must be subcritical.
  if (isnan(d.val)) {
    d = dc;
  }
  d = Length(std::max(dc.val, d.val));
  h = d * a + link->invert(node);
  std::function<bool(Length, Length)> subcritical = [=](Length x, Length h) {
    Length z = S * x + link->invert(node);
    Length d = (h - z) / a;
    if (d > dc) {
      return true;
    } else {
      return false;
    }
  };
  // integrate until the end is reached or the flow becomes supercritical and a
  // hydraulic jump occurs.
  link->HGL = RKF45<Length, Length>(
      std::bind(&PassageLink::hydrualic_slope_subcritical, link,
                std::placeholders::_1, std::placeholders::_2, node),
      h, {Length(0.0), L}, subcritical);
  // get the last entry in the results and return the total energy head at the
  // upstream node.
  h = link->HGL[link->HGL.size() - 1].second;
  Length x = link->HGL[link->HGL.size() - 1].first;
  Length z = link->invert(link->antinode(node)) - S * (L - x);
  Length hv = link->velocity_head(h - z);
  return h + hv;
}

Length frontwater_loss(PassageLink *link, HydraulicNode *node, Length jump_x) {
  // clear recorded HGL profile and initialize data
  link->HGL.clear();
  link->is_supercritical = true;
  if (link->Q.val == 0.0) {
    auto up_node = link->antinode(node);
    node->H = link->invert(node);
    up_node->H = link->invert(up_node);
    link->HGL.push_back({Length(0.0), node->H});
    link->HGL.push_back({Length(link->horizontal_length()), up_node->H});
    return up_node->H;
  }
  Length h = link->get_upstream_water_surface(node);
  Angle S = link->slope(node);
  Length L = link->horizontal_length();
  Dimensionless a = sqrt(S * S + Dimensionless(1.0));
  Length d = (h - link->invert(link->antinode(node))) / a;
  Length dc = link->critical_depth();
  // passage must be steep with starting depth less than or equal to critical
  // depth.
  if (!link->is_steep(node)) {
    link->HGL.push_back({L, Length(nan(""))});
    return Length(nan(""));
  }
  // flow must be supercritical.
  d = Length(std::min(dc.val, d.val));
  h = d * a + link->invert(link->antinode(node));
  std::function<bool(Length, Length)> supercritical = [=](Length x, Length h) {
    Length z = S * x + link->invert(node);
    Length d = (h - z) / a;
    if (d > dc) {
      return false;
    } else {
      return true;
    }
  };
  // integrate until the end is reached or the flow becomes supercritical and a
  // hydraulic jump occurs.
  link->HGL = RKF45<Length, Length>(
      std::bind(&PassageLink::hydrualic_slope_supercritical, link,
                std::placeholders::_1, std::placeholders::_2, node),
      h, {L, jump_x}, supercritical);
  // get the first entry in the results and return the total energy head at the
  // upstream node.
  h = link->HGL[0].second;
  Length z = link->invert(link->antinode(node));
  Length hv = link->velocity_head(h - z);
  return h + hv;
}

Length opening_loss(OpeningLink *link, HydraulicNode *node) {
  // clear recorded HGL profile and initialize data
  link->HGL.clear();
  if (link->Q.val == 0.0) {
    Length H = std::max(node->H, link->elevation);
    link->antinode(node)->H = H;
    link->HGL.push_back({Length(0.0), node->H});
    link->HGL.push_back({Length(0.0), H});
    return H;
  }
  Length h = link->get_water_surface_subcritical(node);
  Length d2 = Length(std::max(h.val - link->elevation.val, 0.0));
  link->HGL.push_back({Length(0.0), h});
  // goal seek Q for d1
  Length dc = critical_depth(link->cross_section_shape.get(), link->Q);
  Length a = std::max(dc, d2);
  Length b = a + Length(0.1);
  Length d1 = find_goal_secant<Length, Flow>(
      abs(link->Q), a, b, 0.00000001,
      std::bind(&OpeningLink::flow, link, std::placeholders::_1, d2));
  h = d1 + link->elevation;
  link->HGL.push_back({Length(0.0), h});
  return h;
}

Length transition_loss(TransitionLink *link, HydraulicNode *node) {
  // clear recorded HGL profile and initialize data
  link->HGL.clear();
  if (link->Q.val == 0.0) {
    auto up_node = link->antinode(node);
    Length H = std::max(node->H, link->invert(up_node));
    up_node->H = H;
    link->HGL.push_back({Length(0.0), node->H});
    link->HGL.push_back({Length(0.0), H});
    return H;
  }
  Length h = link->get_water_surface_subcritical(node);
  Length d2 = h - link->invert(node);
  HydraulicNode *up_node = link->antinode(node);
  Velocity V2 = link->velocity(node, d2);
  Length z2 = link->invert(node);
  Length z1 = link->invert(up_node);

  // objective function to solve for d1.
  std::function<Length(Length)> objective = [=](Length d1) {
    Velocity V1 = link->velocity(up_node, d1);
    if (V1 < V2) {
      // use velocity in the vena contracta
      Dimensionless mu =
          Dimensionless(0.63) +
          Dimensionless(0.37) *
              Dimensionless(pow(V1.val / V2.val, 3.0)); // according to Weisbach
      Velocity V3 = V2 / mu;
      return z2 + d2 + V2 * V2 / (Dimensionless(2.0) * g) +
             link->K * (V3 - V2) * (V3 - V2) / (Dimensionless(2.0) * g) - z1 -
             d1 - V1 * V1 / (Dimensionless(2.0) * g);
    }
    return z2 + d2 + V2 * V2 / (Dimensionless(2.0) * g) +
           link->K * (V1 - V2) * (V1 - V2) / (Dimensionless(2.0) * g) - z1 -
           d1 - V1 * V1 / (Dimensionless(2.0) * g);
  };

  Length d1 = find_goal_secant<Length, Length>(
      Length(0.0), d2, d2 + Length(0.1), 0.0000001, objective);
  // if no solution, set d1 to critical depth
  d1 = isnan(d1.val) || d1 < d2 ? critical_depth(link->shape(up_node), link->Q)
                                : d1;
  // write results to HGL vector
  link->HGL.push_back({Length(0.0), h});
  link->HGL.push_back({Length(0.0), z1 + d1});
  // return total energy head upstream of transition
  return z1 + d1 + link->velocity_head(up_node, d1);
}

// double benching_coefficient(BENCH_CONFIGURATION bench, double D, double d) {
//   double x0 = 1.0, x1 = 2.5, y0, y1;
//   switch (bench) {
//   case BENCH_CONFIGURATION::FLAT:
//     y0 = -0.05;
//     y1 = -0.05;
//     break;
//   case BENCH_CONFIGURATION::DEPRESSED:
//     y0 = 0.00;
//     y1 = 0.00;
//     break;
//   case BENCH_CONFIGURATION::HALF:
//     y0 = -0.85;
//     y1 = -0.05;
//     break;
//   case BENCH_CONFIGURATION::FULL:
//     y0 = -0.93;
//     y1 = -0.25;
//     break;
//   case BENCH_CONFIGURATION::IMPROVED:
//     y0 = -0.98;
//     y1 = -0.60;
//     break;
//   default:
//     y0 = 0.00;
//     y1 = 0.00;
//     break;
//   }
//   if (d <= 0) {
//     d = D;
//   }
//   double x;
//   if (D <= 0.0) {
//     x = 0.0;
//   } else {
//     x = D / d;
//   }
//
//   // If D/d >= 2.5, submerged
//   if (x >= 2.5) {
//     return y1;
//   }
//   // If D/d <= 1.0, unsubmerged
//   if (x <= 1.0) {
//     return y0;
//   }
//   // linearly interpolate y
//   double a = (x - x0) / (x1 - x0);
//   return y0 * (1.0 - a) + y1 * a;
// }
//
// double manhole_loss(double E_ai, double E_i, double Q_0, double D_0,
//                     double invert_0, BENCH_CONFIGURATION bench,
//                     std::vector<double> Q, std::vector<double> invert,
//                     std::vector<double> theta,
//                     std::vector<std::pair<double, double>> &HGL) {
//   // Error handling
//   if (Q_0 <= 0.0 || D_0 <= 0.0) {
//     return E_ai;
//   }
//
//   // find all plunging flows j and non plunging flows k
//   std::vector<double> theta_j{};
//   std::vector<double> Q_j{};
//   std::vector<double> invert_k{};
//   std::vector<double> Q_k{};
//   for (int i = 0; i < invert.size(); i++) {
//     if (invert[i] < E_ai) {
//       // pipe is submerged
//       theta_j.push_back(theta[i]);
//       Q_j.push_back(Q[i]);
//     } else {
//       // pipe is not submerged
//       invert_k.push_back(invert[i]);
//       Q_k.push_back(Q[i]);
//     }
//   }
//
//   // Benching
//   ------------------------------------------------------------------ double
//   C_B = benching_coefficient(bench, D_0, E_ai - invert_0);
//
//   // Angled Inflow
//   ------------------------------------------------------------- double
//   theta_w = 180.0; double Q_j_sum = std::accumulate(std::begin(Q_j),
//   std::end(Q_j), 0.0); if (Q_j_sum > 0.0) {
//     theta_w = std::inner_product(std::begin(Q_j), std::end(Q_j),
//                                  std::begin(theta_j), 0.0) /
//               Q_j_sum;
//   }
//   double C_theta = 4.5 * (Q_j_sum / Q_0) * cos((M_PI / 180.0) * theta_w
//   / 2.0);
//
//   // Plunging
//   ------------------------------------------------------------------
//   std::vector<double> h_k{};
//   for (int k = 0; k < invert_k.size(); k++) {
//     double h = (invert_k[k] - E_ai) / D_0;
//     if (h > 10.0 * D_0) {
//       h = 10.0 * D_0;
//     }
//     h_k.push_back(h);
//   }
//   double C_P =
//       std::inner_product(std::begin(Q_k), std::end(Q_k), std::begin(h_k),
//       0.0) / Q_0;
//
//   // Return results
//   double H_a = std::max(0.0, (C_B + C_theta + C_P) * (E_ai - E_i));
//   HGL.clear();
//   HGL.push_back({0.0, E_ai});
//   HGL.push_back({0.0, E_ai + H_a});
//   return E_ai + H_a;
// }
//
// double bar_screen_loss(HydraulicShape *upstream_channel_shape,
//                        double channel_invert,
//                        HydraulicShape *bar_screen_opening_shape, int n,
//                        double Q, double h) {
//   std::function<double(double)> bar_screen_approach_velocity = [=](double hl)
//   {
//     double upstream_wse = h + hl;
//     double upstream_depth = upstream_wse - channel_invert;
//     double A = upstream_channel_shape->flow_area(upstream_depth);
//     return Q / A;
//   };
//   std::function<double(double)> objective = [=](double hl) {
//     double u = bar_screen_approach_velocity(hl);
//     double d = h - channel_invert; // depth of flow through the bar screen.
//     double q = Q / n;
//     double A = bar_screen_opening_shape->flow_area(d);
//     double v = q / A;
//     return (v * v - u * u) / (0.7 * 2.0 * g) - hl;
//   };
//   double hl = find_goal_secant(0.0, 1.0, 0.5, 0.0000001, objective);
//   double u = bar_screen_approach_velocity(hl);
//   double H =
//       h + hl + (u * u) / (2.0 * g); // upstream energy head in the channel.
//   return H;
// }

} // namespace hazen