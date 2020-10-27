#include "HydraulicLinks.hpp"
#include "HeadLoss.hpp"
#include <algorithm>

namespace hazen {

// Hydraulic Links -------------------------------------------------------------
double HydraulicLink::get_downstream_water_surface() {
  double h = -std::numeric_limits<double>::infinity();
  HydraulicNode *node = Q >= 0.0 ? dn_node : up_node;
  for (auto link : node->links) {
    if ((link->up_node == node && link->Q >= 0.0) ||
        (link->dn_node == node && link->Q < 0.0)) {
      h = std::max(link->HGL[link->HGL.size() - 1].second, h);
    }
  }
  return h == -std::numeric_limits<double>::infinity() ? node->H : h;
}
double HydraulicLink::get_downstream_velocity() {
  HydraulicNode *node = Q >= 0.0 ? up_node : dn_node;
  for (auto &link : node->links) {
    if (auto passage = dynamic_cast<PassageLink *>(link)) {
      double d;
      if (passage->Q >= 0.0) {
        d = passage->HGL[passage->HGL.size() - 1].second - passage->up_inv.z;
      } else {
        d = passage->HGL[passage->HGL.size() - 1].second - passage->dn_inv.z;
      }
      return abs(Q) / passage->cross_section_shape->flow_area(d);
    }
  }
  return 0.0;
}

// Conduits and Channels -------------------------------------------------------
PassageLink::PassageLink(std::unique_ptr<HydraulicShape> cross_section_shape,
                         std::unique_ptr<FrictionMethod> friction_method,
                         vec3 up_inv, vec3 dn_inv, double ds)
    : cross_section_shape(std::move(cross_section_shape)),
      friction_method(std::move(friction_method)), up_inv(up_inv),
      dn_inv(dn_inv), ds(ds) {}
double PassageLink::head_loss() {
  double h = get_downstream_water_surface();
  // if pipe is vertical, use vertical shaft loss.
  if (horizontal_length(dn_inv, up_inv) == 0.0) {
    if (Q >= 0.0) {
      return vertical_shaft_loss(cross_section_shape.get(),
                                 friction_method.get(), up_inv.z, dn_inv.z, Q,
                                 h, HGL);
    } else if (Q < 0.0) {
      return vertical_shaft_loss(cross_section_shape.get(),
                                 friction_method.get(), dn_inv.z, up_inv.z, -Q,
                                 h, HGL);
    }
  }
  double H1 = 0.0;
  double jump_x = nan("");
  std::vector<std::pair<double, double>> passage_backwater_HGL{};
  std::vector<std::pair<double, double>> passage_frontwater_HGL{};
  if (Q >= 0.0) {
    H1 = gvf_backwater_loss(cross_section_shape.get(), friction_method.get(),
                            up_inv, dn_inv, Q, h, ds, jump_x,
                            passage_backwater_HGL);
    if (!isnan(jump_x)) {
      h = get_upstream_water_surface();
      H1 = gvf_frontwater_loss(cross_section_shape.get(), friction_method.get(),
                               up_inv, dn_inv, Q, h, ds, jump_x,
                               passage_frontwater_HGL);
    }
  } else if (Q < 0.0) {
    H1 = gvf_backwater_loss(cross_section_shape.get(), friction_method.get(),
                            dn_inv, up_inv, -Q, h, ds, jump_x,
                            passage_backwater_HGL);
    if (!isnan(jump_x)) {
      h = get_upstream_water_surface();
      H1 = gvf_frontwater_loss(cross_section_shape.get(), friction_method.get(),
                               dn_inv, up_inv, -Q, h, ds, jump_x,
                               passage_frontwater_HGL);
    }
  }
  std::reverse(passage_frontwater_HGL.begin(), passage_frontwater_HGL.end());
  passage_backwater_HGL.insert(passage_backwater_HGL.end(),
                               passage_frontwater_HGL.begin(),
                               passage_frontwater_HGL.end());
  HGL = passage_backwater_HGL;
  return H1;
}
double PassageLink::get_upstream_water_surface() {
  std::vector<double> elevations{};
  HydraulicNode *node = Q >= 0.0 ? up_node : dn_node;
  for (auto link : node->links) {
    if (auto passage = dynamic_cast<PassageLink *>(link)) {
      // regular flow from passage
      if (passage->dn_node == node && passage->Q >= 0.0) {
        double S = slope(passage->dn_inv, passage->up_inv);
        double dn = normal_depth(passage->cross_section_shape.get(),
                                 passage->friction_method.get(), S, passage->Q);
        double dc =
            critical_depth(passage->cross_section_shape.get(), passage->Q);
        bool is_steep = dn < dc ? true : false;
        if (is_steep) {
          double h = passage->get_upstream_water_surface();
          std::vector<std::pair<double, double>> passage_hgl{};
          h = gvf_frontwater_loss(passage->cross_section_shape.get(),
                                  passage->friction_method.get(),
                                  passage->up_inv, passage->dn_inv, passage->Q,
                                  h, ds, 0.0, passage_hgl);
          elevations.push_back(passage_hgl[passage_hgl.size() - 1].second);
        }
      }
      // reverse flow from passage
      if (passage->up_node == node && passage->Q < 0.0) {
        double S = slope(passage->up_inv, passage->dn_inv);
        double dn =
            normal_depth(passage->cross_section_shape.get(),
                         passage->friction_method.get(), S, -passage->Q);
        double dc =
            critical_depth(passage->cross_section_shape.get(), -passage->Q);
        bool is_steep = dn < dc ? true : false;
        if (is_steep) {
          double h = passage->get_upstream_water_surface();
          std::vector<std::pair<double, double>> passage_hgl{};
          h = gvf_frontwater_loss(passage->cross_section_shape.get(),
                                  passage->friction_method.get(),
                                  passage->dn_inv, passage->up_inv, -passage->Q,
                                  h, ds, 0.0, passage_hgl);
          elevations.push_back(passage_hgl[passage_hgl.size() - 1].second);
        }
      }
    }
  }
  if (elevations.size()) {
    return *std::max_element(elevations.begin(), elevations.end());
  }
  double dc = critical_depth(cross_section_shape.get(), abs(Q));
  return Q <= 0.0 ? dn_inv.z + dc : up_inv.z + dc;
}

// Orifices and Weirs ----------------------------------------------------------
OpeningLink::OpeningLink(std::unique_ptr<HydraulicShape> cross_section_shape,
                         double Cd, double elevation, double dy,
                         double percent_open)
    : cross_section_shape(std::move(cross_section_shape)), Cd(Cd),
      elevation(elevation), dy(dy), percent_open(percent_open) {}
double OpeningLink::head_loss() {
  double h = get_downstream_water_surface();
  return opening_loss(cross_section_shape.get(), Cd, elevation, abs(Q), h, dy,
                      HGL, percent_open);
}

// Minor Loss ------------------------------------------------------------------
MinorLink::MinorLink(double K_pos, double K_neg) : K_pos(K_pos), K_neg(K_neg) {}
double MinorLink::head_loss() {
  double V = get_downstream_velocity();
  double h = get_downstream_water_surface();
  double hv = V * V / (2.0 * g);
  double H1 = Q >= 0.0 ? K_pos * hv : K_neg * hv;
  HGL.clear();
  HGL.push_back({0.0, h + hv});
  HGL.push_back({0.0, H1});
  return H1;
}

// Expansion and Contraction ---------------------------------------------------
TransitionLink::TransitionLink(double K) : K(K) {}
double TransitionLink::head_loss() {
  HydraulicShape *up_shape = nullptr, *dn_shape = nullptr;
  double up_inv = -std::numeric_limits<double>::infinity(),
         dn_inv = -std::numeric_limits<double>::infinity();
  if (Q >= 0.0) {
    for (auto &link : up_node->links) {
      if (auto passage = dynamic_cast<PassageLink *>(link)) {
        if (passage->Q >= 0.0) {
          up_inv = passage->dn_inv.z;
        } else {
          up_inv = passage->up_inv.z;
        }
        up_shape = passage->cross_section_shape.get();
      }
    }
    for (auto &link : dn_node->links) {
      if (auto passage = dynamic_cast<PassageLink *>(link)) {
        if (passage->Q >= 0.0) {
          dn_inv = passage->up_inv.z;
        } else {
          dn_inv = passage->dn_inv.z;
        }
        dn_shape = passage->cross_section_shape.get();
      }
    }
  } else {
    for (auto &link : up_node->links) {
      if (auto passage = dynamic_cast<PassageLink *>(link)) {
        if (passage->Q >= 0.0) {
          dn_inv = passage->dn_inv.z;
        } else {
          dn_inv = passage->up_inv.z;
        }
        dn_shape = passage->cross_section_shape.get();
      }
    }
    for (auto &link : dn_node->links) {
      if (auto passage = dynamic_cast<PassageLink *>(link)) {
        if (passage->Q >= 0.0) {
          up_inv = passage->up_inv.z;
        } else {
          up_inv = passage->dn_inv.z;
        }
        up_shape = passage->cross_section_shape.get();
      }
    }
  }
  // if no passage, assume infinite reservoir.
  Rectangle inf_reservoir =
      Rectangle(std::numeric_limits<double>::infinity(),
                std::numeric_limits<double>::infinity(), true);
  if (up_shape == nullptr) {
    up_shape = &inf_reservoir;
  }
  if (dn_shape == nullptr) {
    dn_shape = &inf_reservoir;
  }
  double h = get_downstream_water_surface();
  double H1 =
      transition_loss(up_shape, dn_shape, up_inv, dn_inv, abs(Q), K, h, HGL);
  return H1;
}

// Manhole losses --------------------------------------------------------------
ManholeLink::ManholeLink(double elevation, BENCH_CONFIGURATION bench_config)
    : elevation(elevation), bench_config(bench_config) {}
double ManholeLink::head_loss() {
  // find flow direction of all connecting pipes and bind them to the correct
  // nodes

  // find the parameters for the head loss equations

  // clac new H upstream.
}

double NullLink::head_loss() { return 0.0; }

} // namespace hazen