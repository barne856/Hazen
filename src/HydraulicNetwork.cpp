#include "HydraulicNetwork.hpp"
#include <cmath>

namespace hazen {

HydraulicNode::HydraulicNode() : E(nan("")) {}
double HydraulicNode::continuity() {
  double Q = 0.;
  for (const auto &link : links) {
    if (link->dn_node == this) {
      Q += link->Q;
    }
    if (link->up_node == this) {
      Q -= link->Q;
    }
  }
  return Q;
}

HydraulicLink::HydraulicLink()
    : Q(nan("")), up_node(nullptr), dn_node(nullptr) {}

void HydraulicNetwork::add_component(
    std::shared_ptr<HydraulicComponent> component) {
  components.insert(component);
}
} // namespace hazen