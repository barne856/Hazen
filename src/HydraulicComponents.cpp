#include "HydraulicComponents.hpp"

namespace hazen {
// Outfall ---------------------------------------------------------------------
Outfall::Outfall() : outfall_node(std::make_shared<HydraulicNode>()) {}
void Outfall::bind(std::pair<HydraulicComponent *, unsigned int> bind_point,
                   unsigned int binding_index) {
  std::shared_ptr<HydraulicLink> link;
  if (dynamic_cast<Passage *>(bind_point.first)) {
    link = std::make_shared<TransitionLink>(1.0);
  } else {
    link = std::make_shared<NullLink>();
  }
  bind_components({this, binding_index}, bind_point, link);
}
HydraulicNode *Outfall::get_binding_node(unsigned int binding_index) {
  return outfall_node.get();
}

// Manhole ---------------------------------------------------------------------
Manhole::Manhole(double invert, BENCH_CONFIGURATION bench_config)
    : invert(invert), bench_config(bench_config),
      manhole_link(std::make_shared<ManholeLink>()) {}

// Storage ---------------------------------------------------------------------
Storage::Storage(std::vector<std::pair<double, double>> storage_curve,
                 double elevation)
    : storage_curve(storage_curve), elevation(elevation),
      storage_node(std::make_shared<HydraulicNode>()) {}
void Storage::bind(std::pair<HydraulicComponent *, unsigned int> bind_point,
                   unsigned int binding_index) {
  std::shared_ptr<HydraulicLink> link;
  if (dynamic_cast<Passage *>(bind_point.first)) {
    link = std::make_shared<TransitionLink>(1.0);
  } else {
    link = std::make_shared<NullLink>();
  }
  bind_components({this, binding_index}, bind_point, link);
}
HydraulicNode *Storage::get_binding_node(unsigned int binding_index) {
  return storage_node.get();
}
} // namespace hazen
