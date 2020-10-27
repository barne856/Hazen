#include "HydraulicNetwork.hpp"
#include <algorithm>

namespace hazen {

HydraulicNode::HydraulicNode() : H(0.0) {}
HydraulicNode::~HydraulicNode() {
  for (auto &link : links) {
    if (link->dn_node == this) {
      link->dn_node = nullptr;
    }
    if (link->up_node == this) {
      link->up_node = nullptr;
    }
  }
}
double HydraulicNode::continuity() {
  double Q = 0.0;
  for (const auto &link : links) {
    if (link->dn_node == this) {
      Q += link->Q;
    }
    if (link->up_node == this) {
      Q -= link->Q;
    }
  }
  for (const auto &flow : point_flows) {
    Q += flow;
  }
  return Q;
}

HydraulicLink::HydraulicLink() : Q(0.0), up_node(nullptr), dn_node(nullptr) {}
HydraulicLink::~HydraulicLink() {
  set_up_node(nullptr);
  set_dn_node(nullptr);
}
void HydraulicLink::set_up_node(HydraulicNode *node) {
  if (up_node) {
    up_node->links.erase(this);
  }
  up_node = node;
  if (up_node) {
    up_node->links.insert(this);
  }
}
void HydraulicLink::set_dn_node(HydraulicNode *node) {
  if (dn_node) {
    dn_node->links.erase(this);
  }
  dn_node = node;
  if (dn_node) {
    dn_node->links.insert(this);
  }
}

HydraulicComponent::~HydraulicComponent() {
  for (const auto &binding_point : binding_points) {
    unbind(binding_point.first);
  }
}
void HydraulicComponent::unbind(unsigned int binding_index) {
  // remove binding links and binding points from the bound component
  binding_points[binding_index].first->binding_links.erase(
      binding_points[binding_index].second);
  binding_points[binding_index].first->binding_points.erase(
      binding_points[binding_index].second);
  // remove binding links and binding points from this component
  binding_links.erase(binding_index);
  binding_points.erase(binding_index);
}
void HydraulicComponent::bind_components(
    std::pair<HydraulicComponent *, unsigned int> bind_point_1,
    std::pair<HydraulicComponent *, unsigned int> bind_point_2,
    std::shared_ptr<HydraulicLink> link) {
  link->set_up_node(bind_point_1.first->get_binding_node(bind_point_1.second));
  link->set_dn_node(bind_point_2.first->get_binding_node(bind_point_2.second));
  bind_point_1.first->binding_links.insert({bind_point_1.second, link});
  bind_point_1.first->binding_points.insert(
      {bind_point_1.second, bind_point_2});
  bind_point_2.first->binding_links.insert({bind_point_2.second, link});
  bind_point_2.first->binding_points.insert(
      {bind_point_2.second, bind_point_1});
}

} // namespace hazen