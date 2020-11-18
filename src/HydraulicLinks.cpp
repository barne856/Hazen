#include "HydraulicLinks.hpp"
#include "HeadLoss.hpp"
#include "HydraulicUtil.hpp"
#include <algorithm>

namespace hazen {

// Conduits and Channels -------------------------------------------------------
PassageLink::PassageLink(std::shared_ptr<HydraulicShape> cross_section_shape,
                         std::shared_ptr<FrictionMethod> friction_method,
                         std::pair<vec3, vec3> alignment)
    : cross_section_shape(cross_section_shape),
      friction_method(friction_method), alignment(alignment) {}
double PassageLink::velocity_head(double depth) {
  double A = cross_section_shape->flow_area(depth);
  double V = isinf(A) || isnan(A) || A == 0.0 ? 0.0 : abs(Q) / A;
  return V * V / (2.0 * g);
}
double PassageLink::invert(HydraulicNode *node) {
  if (this->node<0>() == node) {
    return alignment.first.z;
  }
  if (this->node<1>() == node) {
    return alignment.second.z;
  }
  return nan("");
}
double PassageLink::length() {
  return hazen::length(alignment.first, alignment.second);
}
double PassageLink::horizontal_length() {
  return hazen::horizontal_length(alignment.first, alignment.second);
}
double PassageLink::slope(HydraulicNode *node) {
  vec3 up = {nan(""), nan(""), nan("")}, dn = {nan(""), nan(""), nan("")};
  if (this->node<0>() == node) {
    up = alignment.second;
    dn = alignment.first;
  } else if (this->node<1>() == node) {
    up = alignment.first;
    dn = alignment.second;
  }
  return hazen::slope(dn, up);
}
double PassageLink::critical_depth() {
  return hazen::critical_depth(cross_section_shape.get(), abs(Q));
}
double PassageLink::normal_depth(HydraulicNode *node) {
  return hazen::normal_depth(cross_section_shape.get(), friction_method.get(),
                             slope(node), abs(Q));
}
bool PassageLink::is_steep(HydraulicNode *node) {
  return normal_depth(node) < critical_depth() ? true : false;
}
double PassageLink::friction_slope(double depth) {
  return friction_method->friction_slope(cross_section_shape.get(), abs(Q),
                                         depth);
}
double PassageLink::hydrualic_slope_subcritical(double x, double h,
                                                HydraulicNode *node) {
  double S = slope(node);
  double a = sqrt(S * S + 1); // large angle correction.
  double d = (h - (invert(node) + S * x)) / a;
  double dc = critical_depth();
  double TOL = 0.01; // tolerance for numerical stability
  if (d < dc) {
    // force subcritical flow
    d = dc + TOL;
  }
  double Sf = friction_slope(d);
  double Fr = cross_section_shape->froude(Q, d);
  double Fr_squared = Fr * Fr;
  if (1.0 - Fr_squared < TOL) {
    // force subcritical flow
    Fr_squared = 1.0 - TOL;
  }
  return (Sf * a - S) / (1.0 - Fr_squared) + S;
}
double PassageLink::hydrualic_slope_supercritical(double x, double h,
                                                  HydraulicNode *node) {
  double S = slope(node);
  double a = sqrt(S * S + 1); // large angle correction.
  double d = (h - (invert(node) + S * x)) / a;
  double dc = critical_depth();
  double TOL = 0.01; // tolerance for numerical stability
  if (d < TOL) {
    return S;
  }
  if (d > dc) {
    // force supercritical flow
    d = dc - TOL;
  }
  double Sf = friction_slope(d);
  double Fr = cross_section_shape->froude(Q, d);
  double Fr_squared = Fr * Fr;
  if (1.0 - Fr_squared > -TOL) {
    // force supercritical flow
    Fr_squared = 1.0 + TOL;
  }
  return (Sf * a - S) / (1.0 - Fr_squared) + S;
}

double PassageLink::head_loss(HydraulicNode *node) {
  HydraulicNode *up_node = antinode(node);
  // vertical pipes
  if (horizontal_length() == 0.0) {
    if (invert(node) <= invert(up_node)) {
      double H = vertical_drop_loss(this, node);
      up_node->H = H;
      return H;
    } else {
      double H = vertical_rise_loss(this, node);
      up_node->H = H;
      return H;
    }
  }
  // test for hydraulic jump with backwater calcs
  double H = backwater_loss(this, node);
  if (isnan(H)) {
    auto temp_HGL = HGL;
    double jump_x = temp_HGL[temp_HGL.size() - 1].first;
    temp_HGL.pop_back();
    H = frontwater_loss(this, node, jump_x);
    for (auto entry : temp_HGL) {
      HGL.push_back(entry);
    }
    std::sort(HGL.begin(), HGL.end());
  }
  up_node->H = H;
  return H;
}
double PassageLink::get_water_surface_subcritical(HydraulicNode *node) {
  std::function<double(double)> objective = [=](double h) {
    double S = slope(node);
    double a = sqrt(S * S + 1.0);
    double d = (h - invert(node)) / a;
    return h + velocity_head(d) - node->H;
  };
  double h =
      find_goal_secant(0.0, node->H, node->H - 0.001, 0.000001, objective);
  // nan is returned to signal an H that is too small for subcritical flow.
  return h;
}
double PassageLink::get_water_surface_supercritical(HydraulicNode *node) {
  std::function<double(double)> objective = [=](double h) {
    double S = slope(node);
    double a = sqrt(S * S + 1.0);
    double d = (h - invert(node)) / a;
    return h + velocity_head(d) - node->H;
  };
  double h = find_goal_secant(0.0, invert(node) + 0.001, invert(node) + 0.0015,
                              0.000001, objective);
  // nan is returned to signal an H that is too small for supercritical flow.
  return h;
}
double PassageLink::get_upstream_water_surface(HydraulicNode *node) {
  HydraulicNode *up_node = antinode(node);
  std::vector<double> elevations{};
  // loop through all passages steeply flowing into the upstream node.
  for (auto i : up_node->in_links()) {
    if (auto passage =
            std::dynamic_pointer_cast<PassageLink>(up_node->links[i])) {
      if (passage->is_steep(up_node)) {
        // record water surface in a vector
        frontwater_loss(passage.get(), up_node);
        elevations.push_back(passage->HGL[passage->HGL.size() - 1].second);
      }
    }
  }
  // if there are steep inflowing links, return the max h recorded.
  if (elevations.size()) {
    return *std::max_element(elevations.begin(), elevations.end());
  }
  // else return the critical depth upstream.
  return critical_depth() + invert(up_node);
}

// Orifices and Weirs ----------------------------------------------------------
OpeningLink::OpeningLink(std::shared_ptr<HydraulicShape> cross_section_shape,
                         double Cd, double elevation, double percent_open)
    : cross_section_shape(std::move(cross_section_shape)), Cd(Cd),
      elevation(elevation), percent_open(percent_open) {}
double OpeningLink::velocity(double y, double d1, double d2) {
  double D = cross_section_shape->get_max_depth() * percent_open;
  if (0.0 <= y && y < std::max(0.0, d2)) {
    return sqrt(2.0 * g * (d1 - d2));
  } else if (std::max(0.0, d2) <= y && y <= std::min(d1, D)) {
    return sqrt(2.0 * g * (d1 - y));
  } else {
    return 0.0;
  }
}
double OpeningLink::flow_step(double y, double Q, double d1, double d2) {
  return Cd * velocity(y, d1, d2) * cross_section_shape->top_width(y);
}
double OpeningLink::flow(double d1, double d2) {
  double D = cross_section_shape->get_max_depth() * percent_open;
  auto result =
      RKF45(std::bind(&OpeningLink::flow_step, this, std::placeholders::_1,
                      std::placeholders::_2, d1, d2),
            0.0, {0.0, D});
  return result.back().second;
}
double OpeningLink::get_water_surface_subcritical(HydraulicNode *node) {
  return node->H;
}
double OpeningLink::head_loss(HydraulicNode *node) {
  HydraulicNode *up_node = antinode(node);
  double H = opening_loss(this, node);
  up_node->H = H;
  return H;
}

// Minor Loss ------------------------------------------------------------------
MinorLink::MinorLink(double K_pos, double K_neg,
                     std::shared_ptr<HydraulicShape> cross_section_shape,
                     double invert)
    : K_pos(K_pos), K_neg(K_neg), cross_section_shape(cross_section_shape),
      invert(invert) {}
double MinorLink::velocity_head(double depth) {
  double A = cross_section_shape->flow_area(depth);
  double V = isinf(A) || isnan(A) || A == 0.0 ? 0.0 : abs(Q) / A;
  return V * V / (2.0 * g);
}
double MinorLink::get_water_surface_subcritical(HydraulicNode *node) {
  std::function<double(double)> objective = [=](double h) {
    double d = h - invert;
    return h + velocity_head(d) - node->H;
  };
  double h =
      find_goal_secant(0.0, node->H, node->H - 0.001, 0.000001, objective);
  // check if nan is returned. If so, return critical depth elevation.
  if (isnan(h)) {
    return invert + critical_depth(cross_section_shape.get(), Q);
  }
  return h;
}
double MinorLink::head_loss(HydraulicNode *node) {
  double h = get_water_surface_subcritical(node);
  double hv = velocity_head(h - invert);
  HydraulicNode *up_node = antinode(node);
  double K = Q >= 0 ? K_pos : K_neg;
  double hl = K * hv;
  double H = node->H + hl;
  up_node->H = H;
  return H;
}

// Expansion and Contraction ---------------------------------------------------
TransitionLink::TransitionLink(
    double K,
    std::pair<std::shared_ptr<HydraulicShape>, std::shared_ptr<HydraulicShape>>
        shapes,
    std::pair<double, double> inverts)
    : K(K), shapes(shapes), inverts(inverts) {}
double TransitionLink::get_water_surface_subcritical(HydraulicNode *node) {
  std::function<double(double)> objective = [=](double h) {
    double d = h - invert(node);
    return h + velocity_head(node, d) - node->H;
  };
  double h =
      find_goal_secant(0.0, node->H, node->H - 0.001, 0.000001, objective);
  // check if nan is returned. If so, return critical depth elevation.
  if (isnan(h)) {
    return invert(node) + critical_depth(shape(node), Q);
  }
  return h;
}
double TransitionLink::invert(HydraulicNode *node) {
  if (this->node<0>() == node) {
    return inverts.first;
  }
  if (this->node<1>() == node) {
    return inverts.second;
  }
  return nan("");
}
HydraulicShape *TransitionLink::shape(HydraulicNode *node) {
  if (this->node<0>() == node) {
    return shapes.first.get();
  }
  if (this->node<1>() == node) {
    return shapes.second.get();
  }
  return nullptr;
}
double TransitionLink::velocity_head(HydraulicNode *node, double depth) {
  double V = velocity(node, depth);
  return V * V / (2.0 * g);
}
double TransitionLink::velocity(HydraulicNode *node, double depth) {
  HydraulicShape *the_shape = shape(node);
  if (the_shape) {
    double A = shape(node)->flow_area(depth);
    return isinf(A) || isnan(A) || A == 0.0 ? 0.0 : abs(Q) / A;
  }
  return 0.0;
}
double TransitionLink::head_loss(HydraulicNode *node) {
  HydraulicNode *up_node = antinode(node);
  double H = transition_loss(this, node);
  up_node->H = H;
  return H;
}

} // namespace hazen