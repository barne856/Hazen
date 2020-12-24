#include "Hazen/HydraulicLinks.hpp"
#include "Hazen/HeadLoss.hpp"
#include <algorithm>
#include <iostream>

namespace hazen {

// Conduits and Channels -------------------------------------------------------
PassageLink::PassageLink(std::shared_ptr<HydraulicShape> cross_section_shape,
                         std::shared_ptr<FrictionMethod> friction_method,
                         std::pair<Vec<Length>, Vec<Length>> alignment)
    : cross_section_shape(cross_section_shape),
      friction_method(friction_method), alignment(alignment) {}
Length PassageLink::velocity_head(Length depth) const {
  Area A = cross_section_shape->flow_area(depth);
  Velocity V =
      isinf(A.val) || isnan(A.val) || A.val == 0.0 ? Velocity(0.0) : abs(Q) / A;
  return (V * V) / (Dimensionless(2.0) * g);
}
Length PassageLink::invert(HydraulicNode *node) const {
  if (this->node(0) == node) {
    return alignment.first(2);
  }
  if (this->node(1) == node) {
    return alignment.second(2);
  }
  return Length(nan(""));
}
Length PassageLink::length() const {
  return norm(alignment.first - alignment.second);
}
Length PassageLink::horizontal_length() const {
  return hazen::horizontal_length(alignment.first, alignment.second);
}
Angle PassageLink::slope(HydraulicNode *node) const {
  if (this->node(0) == node) {
    return hazen::slope(alignment.first, alignment.second);
  } else {
    return hazen::slope(alignment.second, alignment.first);
  }
}
Length PassageLink::critical_depth() const {
  return hazen::critical_depth(cross_section_shape.get(), abs(Q));
}
Length PassageLink::normal_depth(HydraulicNode *node) const {
  return hazen::normal_depth(cross_section_shape.get(), friction_method.get(),
                             this->slope(node), abs(Q));
}
bool PassageLink::is_steep(HydraulicNode *node) const {
  return normal_depth(node) < critical_depth() ? true : false;
}
Angle PassageLink::friction_slope(Length depth) const {
  return friction_method->friction_slope(cross_section_shape.get(), abs(Q),
                                         depth);
}
Angle PassageLink::hydrualic_slope_subcritical(Length x, Length h,
                                               HydraulicNode *node) const {
  Angle S = slope(node);
  Dimensionless a = sqrt(S * S + Dimensionless(1.0)); // large angle correction.
  Length d = (h - (invert(node) + S * x)) / a;
  Length dc = critical_depth();
  Dimensionless TOL = Dimensionless(0.01); // tolerance for numerical stability
  if (d < dc) {
    // force subcritical flow
    d = dc + Length(TOL.val);
  }
  Angle Sf = friction_slope(d);
  Dimensionless Fr = cross_section_shape->froude(Q, d);
  Dimensionless Fr_squared = Fr * Fr;
  if (Dimensionless(1.0) - Fr_squared < TOL) {
    // force subcritical flow
    Fr_squared = Dimensionless(1.0) - TOL;
  }
  return (Sf * a - S) / (Dimensionless(1.0) - Fr_squared) + S;
}
Angle PassageLink::hydrualic_slope_supercritical(Length x, Length h,
                                                 HydraulicNode *node) const {
  Angle S = slope(node);
  Dimensionless a = sqrt(S * S + Dimensionless(1.0)); // large angle correction.
  Length d = (h - (invert(node) + S * x)) / a;
  Length dc = critical_depth();
  Dimensionless TOL = Dimensionless(0.01); // tolerance for numerical stability
  if (d.val < TOL.val) {
    return S;
  }
  if (d > dc) {
    // force supercritical flow
    d = dc - Length(TOL.val);
  }
  Angle Sf = friction_slope(d);
  Dimensionless Fr = cross_section_shape->froude(Q, d);
  Dimensionless Fr_squared = Fr * Fr;
  if ((Dimensionless(1.0) - Fr_squared).val > -TOL.val) {
    // force supercritical flow
    Fr_squared = Dimensionless(1.0) + TOL;
  }
  return (Sf * a - S) / (Dimensionless(1.0) - Fr_squared) + S;
}

Length PassageLink::head_loss(HydraulicNode *node) {
  HydraulicNode *up_node = antinode(node);
  // vertical pipes
  if (horizontal_length() == Length(0.0)) {
    if (invert(node) <= invert(up_node)) {
      Length H = vertical_drop_loss(this, node);
      up_node->H = H;
      return H;
    } else {
      Length H = vertical_rise_loss(this, node);
      up_node->H = H;
      return H;
    }
  }
  // test for hydraulic jump with backwater calcs
  Length H = backwater_loss(this, node);
  if (isnan(H.val)) {
    auto temp_HGL = HGL;
    Length jump_x = temp_HGL[temp_HGL.size() - 1].first;
    temp_HGL.pop_back();
    H = frontwater_loss(this, node, jump_x);
    for (auto entry : temp_HGL) {
      HGL.push_back(entry);
    }
    std::sort(HGL.begin(), HGL.end());
  }
  up_node->H = H;
  // if (isnan(H.val)) {
  //   std::cout << "Whoopsie" << std::endl;
  // }
  return H;
}
Length PassageLink::get_water_surface_subcritical(HydraulicNode *node) const {
  std::function<Length(Length)> objective = [=](Length h) {
    Angle S = slope(node);
    Dimensionless a = sqrt(S * S + Dimensionless(1.0));
    Length d = (h - invert(node)) / a;
    return h + velocity_head(d) - node->H;
  };
  Length h = find_goal_secant<Length, Length>(
      Length(0.0), node->H, node->H - Length(0.001), 0.000001, objective);
  // nan is returned to signal an H that is too small for subcritical flow.
  return h;
}
Length PassageLink::get_water_surface_supercritical(HydraulicNode *node) const {
  std::function<Length(Length)> objective = [=](Length h) {
    Angle S = slope(node);
    Dimensionless a = sqrt(S * S + Dimensionless(1.0));
    Length d = (h - invert(node)) / a;
    return h + velocity_head(d) - node->H;
  };
  Length h = find_goal_secant<Length, Length>(
      Length(0.0), invert(node) + Length(0.001), invert(node) + Length(0.0015),
      0.000001, objective);
  // nan is returned to signal an H that is too small for supercritical flow.
  return h;
}
Length PassageLink::get_upstream_water_surface(HydraulicNode *node) const {
  HydraulicNode *up_node = antinode(node);
  std::vector<Length> elevations{};
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
                         Dimensionless Cd, Length elevation,
                         Dimensionless percent_open)
    : cross_section_shape(std::move(cross_section_shape)), Cd(Cd),
      elevation(elevation), percent_open(percent_open) {}
Velocity OpeningLink::velocity(Length y, Length d1, Length d2) const {
  Length D = cross_section_shape->get_max_depth() * percent_open;
  if (0.0 <= y.val && y.val < std::max(0.0, d2.val)) {
    return sqrt(Dimensionless(2.0) * g * (d1 - d2));
  } else if (std::max(0.0, d2.val) <= y.val &&
             y.val <= std::min(d1.val, D.val)) {
    return sqrt(Dimensionless(2.0) * g * (d1 - y));
  } else {
    return Velocity(0.0);
  }
}
quotient_type<Area, Time> OpeningLink::flow_step(Length y, Flow Q, Length d1,
                                                 Length d2) const {
  return Cd * velocity(y, d1, d2) * cross_section_shape->top_width(y);
}
Flow OpeningLink::flow(Length d1, Length d2) const {
  Length D = cross_section_shape->get_max_depth() * percent_open;
  auto result = RKF45<Length, Flow>(std::bind(&OpeningLink::flow_step, this,
                                              std::placeholders::_1,
                                              std::placeholders::_2, d1, d2),
                                    Flow(0.0), {Length(0.0), D});
  return result.back().second;
}
Length OpeningLink::get_water_surface_subcritical(HydraulicNode *node) const {
  return node->H;
}
Length OpeningLink::head_loss(HydraulicNode *node) {
  HydraulicNode *up_node = antinode(node);
  Length H = opening_loss(this, node);
  up_node->H = H;
  return H;
}

// Minor Loss ------------------------------------------------------------------
MinorLink::MinorLink(Dimensionless K_pos, Dimensionless K_neg,
                     std::shared_ptr<HydraulicShape> cross_section_shape,
                     Length invert)
    : K_pos(K_pos), K_neg(K_neg), cross_section_shape(cross_section_shape),
      invert(invert) {}
Length MinorLink::velocity_head(Length depth) const {
  Area A = cross_section_shape->flow_area(depth);
  Velocity V =
      isinf(A.val) || isnan(A.val) || A.val == 0.0 ? Velocity(0.0) : abs(Q) / A;
  return V * V / (Dimensionless(2.0) * g);
}
Length MinorLink::get_water_surface_subcritical(HydraulicNode *node) const {
  std::function<Length(Length)> objective = [=](Length h) {
    Length d = h - invert;
    return h + velocity_head(d) - node->H;
  };
  Length h = find_goal_secant<Length, Length>(
      Length(0.0), node->H, node->H - Length(0.001), 0.000001, objective);
  // check if nan is returned. If so, return critical depth elevation.
  if (isnan(h.val)) {
    return invert + critical_depth(cross_section_shape.get(), Q);
  }
  return h;
}
Length MinorLink::head_loss(HydraulicNode *node) {
  Length h = get_water_surface_subcritical(node);
  Length hv = velocity_head(h - invert);
  HydraulicNode *up_node = antinode(node);
  Dimensionless K = Q.val >= 0.0 ? K_pos : K_neg;
  Length hl = K * hv;
  Length H = node->H + hl;
  up_node->H = H;
  return H;
}

// Expansion and Contraction ---------------------------------------------------
TransitionLink::TransitionLink(
    Dimensionless K,
    std::pair<std::shared_ptr<HydraulicShape>, std::shared_ptr<HydraulicShape>>
        shapes,
    std::pair<Length, Length> inverts)
    : K(K), shapes(shapes), inverts(inverts) {}
Length
TransitionLink::get_water_surface_subcritical(HydraulicNode *node) const {
  std::function<Length(Length)> objective = [=](Length h) {
    Length d = h - invert(node);
    return h + velocity_head(node, d) - node->H;
  };
  Length h = find_goal_secant<Length, Length>(
      Length(0.0), node->H, node->H - Length(0.001), 0.000001, objective);
  // check if nan is returned. If so, return critical depth elevation.
  if (isnan(h.val)) {
    return invert(node) + critical_depth(shape(node), Q);
  }
  return h;
}
Length TransitionLink::invert(HydraulicNode *node) const {
  if (this->node(0) == node) {
    return inverts.first;
  }
  if (this->node(1) == node) {
    return inverts.second;
  }
  return Length(nan(""));
}
HydraulicShape *TransitionLink::shape(HydraulicNode *node) const {
  if (this->node(0) == node) {
    return shapes.first.get();
  }
  if (this->node(1) == node) {
    return shapes.second.get();
  }
  return nullptr;
}
Length TransitionLink::velocity_head(HydraulicNode *node, Length depth) const {
  Velocity V = velocity(node, depth);
  return V * V / (Dimensionless(2.0) * g);
}
Velocity TransitionLink::velocity(HydraulicNode *node, Length depth) const {
  HydraulicShape *the_shape = shape(node);
  if (the_shape) {
    Area A = shape(node)->flow_area(depth);
    return isinf(A.val) || isnan(A.val) || A.val == 0.0 ? Velocity(0.0)
                                                        : abs(Q) / A;
  }
  return Velocity(0.0);
}
Length TransitionLink::head_loss(HydraulicNode *node) {
  HydraulicNode *up_node = antinode(node);
  Length H = transition_loss(this, node);
  up_node->H = H;
  return H;
}

} // namespace hazen