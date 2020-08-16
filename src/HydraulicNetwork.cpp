#include "HydraulicNetwork.hpp"

namespace hazen {
void HydraulicNetwork::add_point_flow(std::shared_ptr<PointFlow> point_flow) {
  point_flows.insert(point_flow);
}
void HydraulicNetwork::add_varied_flow(
    std::shared_ptr<VariedFlow> varied_flow) {
  varied_flows.insert(varied_flow);
}
void HydraulicNetwork::add_component(
    std::shared_ptr<HydraulicComponent> component) {
  components.insert(component);
}
} // namespace hazen