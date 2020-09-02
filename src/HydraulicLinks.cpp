#include "HydraulicLinks.hpp"
#include <cmath>

namespace hazen {
MinorLossLink::MinorLossLink() : K_pos(nan("")), K_neg(nan("")) {}
PassageLink::PassageLink()
    : cross_section_shape(nullptr), friction_method(nullptr) {}
OpeningLink::OpeningLink()
    : Cd(nan("")), elevation(nan("")), cross_section_shape(nullptr) {}
void OpeningLink::set_cross_section(std::shared_ptr<HydraulicShape> shape) {
  cross_section_shape = shape;
}
PumpLink::PumpLink() : elevation(nan("")) {}
SlideGateLink::SlideGateLink()
    : elevation(nan("")), percent_open(1.0), is_closed(false) {}
StopValveLink::StopValveLink() : elevation(nan("")), is_closed(false) {}
CheckValveLink::CheckValveLink() : elevation(nan("")) {}
BaffleLink::BaffleLink() {}
} // namespace hazen