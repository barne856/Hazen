#include "HydraulicLinks.hpp"
#include "HeadLoss.hpp"
#include <cmath>

namespace hazen {

MinorLink::MinorLink() : K(nan("")), downstream_passage(nullptr) {}
double MinorLink::head_loss() {
  if (downstream_passage) {
    double V = downstream_passage->get_up_node_velocity();
    return K * V * V / (2 * g);
  }
  return 0.0;
}
TransitionLink::TransitionLink()
    : K(nan("")), downstream_passage(nullptr), upstream_passage(nullptr) {}
double TransitionLink::head_loss() {
  if (downstream_passage && upstream_passage) {
    double V_dn = downstream_passage->get_up_node_velocity();
    double V_up = upstream_passage->get_dn_node_velocity();
    if (V_dn < V_up) {
      // expansion loss assumed
      std::function<double(double)> F = [=](double h_up) -> double {
        up_node->H = h_up;
        double V = upstream_passage->get_dn_node_velocity();
        return up_node->H - dn_node->H + (V * V - V_dn * V_dn) / (2.0 * g) -
               K * (V * V_dn) * (V * V_dn) / (2.0 * g);
      };
      up_node->H =
          find_goal_secant(0.0, dn_node->H, dn_node->H + 0.1, 0.00001, F);
      V_dn = downstream_passage->get_up_node_velocity();
      V_up = upstream_passage->get_dn_node_velocity();
      double dh = up_node->H - dn_node->H;
      return dh + (V_up * V_up - V_dn * V_dn) / (2.0 * g);
    } else if (V_dn > V_up) {
      // contraction loss assumed
      std::function<double(double)> F = [=](double h_up) -> double {
        up_node->H = h_up;
        double V = upstream_passage->get_dn_node_velocity();
        double A_dn = downstream_passage->get_up_node_area();
        double A_up = upstream_passage->get_dn_node_area();
        double mu = 0.63 + 0.37 * pow(A_dn / A_up, 3.0);
        return up_node->H - dn_node->H + (V * V - V_dn * V_dn) / (2.0 * g) -
               K * (1.0 / 2.0 / g) * pow(1.0 / mu - 1.0, 2.0) *
                   pow(A_up / A_dn, 2.0) * V_up * V_up;
      };
      up_node->H =
          find_goal_secant(0.0, dn_node->H, dn_node->H + 0.1, 0.00001, F);
      V_dn = downstream_passage->get_up_node_velocity();
      V_up = upstream_passage->get_dn_node_velocity();
      double dh = up_node->H - dn_node->H;
      return dh + (V_up * V_up - V_dn * V_dn) / (2.0 * g);
    }
    return 0.0;
  }
  return 0.0;
}
PassageLink::PassageLink()
    : cross_section_shape(nullptr), friction_method(nullptr) {}
void PassageLink::set_cross_section(std::shared_ptr<HydraulicShape> shape) {
  cross_section_shape = shape;
}
void PassageLink::set_friction_method(
    std::shared_ptr<FrictionMethod> friction_method) {
  this->friction_method = friction_method;
}
double PassageLink::get_up_node_velocity() { return Q / get_up_node_area(); }
double PassageLink::get_dn_node_velocity() { return Q / get_dn_node_area(); }
double PassageLink::get_up_node_area() {
  if (up_node) {
    // most upstream point of passage
    vec3 p1 = invert_alignment[invert_alignment.size() - 1];
    // second most upstream point of passage
    vec3 p2 = invert_alignment[invert_alignment.size() - 2];
    // vertical depth to invert at upstream end of passage
    double y = up_node->H - p1.z;
    // Slope of upstream end of passage
    double S = slope(p2, p1);
    // flow depth in passage
    double d;
    if (S == std::numeric_limits<double>::infinity() ||
        S == -std::numeric_limits<double>::infinity()) {
      d = cross_section_shape->get_max_depth();
    } else {
      d = y / sqrt(S * S + 1);
    }
    // velocity at upstream end of passage
    return cross_section_shape->flow_area(d);
  }
  return nan("");
}
double PassageLink::get_dn_node_area() {
  if (dn_node) {
    // second most downstream point of passage
    vec3 p1 = invert_alignment[1];
    // most downstream point of passage
    vec3 p2 = invert_alignment[0];
    // vertical depth to invert at downstream end of passage
    double y = up_node->H - p1.z;
    // Slope of udownstreampstream end of passage
    double S = slope(p2, p1);
    // flow depth in passage
    double d;
    if (S == std::numeric_limits<double>::infinity() ||
        S == -std::numeric_limits<double>::infinity()) {
      d = cross_section_shape->get_max_depth();
    } else {
      d = y / sqrt(S * S + 1);
    }
    // velocity at downstream end of passage
    return cross_section_shape->flow_area(d);
  }
  return nan("");
}
double PassageLink::head_loss() {
  double jump_x = nan("");
  double h_new = dn_node->H;
  for (int i = 0; i < invert_alignment.size() - 1; i++) {
    vec3 dn_inv = invert_alignment[i];
    vec3 up_inv = invert_alignment[i + 1];
    double loss =
        gvf_backwater_loss(cross_section_shape.get(), friction_method.get(),
                           up_inv, dn_inv, Q, h_new, 0.001, jump_x);
    h_new += loss;
    if (!isnan(jump_x)) {
      double S = slope(dn_inv, up_inv);
      double a = sqrt(S * S + 1);
      h_new = up_inv.z + a * critical_depth(cross_section_shape.get(), Q);
    }
  }
  return h_new - dn_node->H;
}
OpeningLink::OpeningLink()
    : Cd(nan("")), elevation(nan("")), cross_section_shape(nullptr) {}
double OpeningLink::head_loss() {
  return opening_loss(cross_section_shape.get(), Cd, elevation, Q, dn_node->H,
                      0.0001);
}
void OpeningLink::set_cross_section(std::shared_ptr<HydraulicShape> shape) {
  cross_section_shape = shape;
}
} // namespace hazen