#ifndef HYDRAULICNETWORK
#define HYDRAULICNETWORK

// Standard Library
#include <memory>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace hazen {
// Forward declare HydraulicLink for use as a pointer in HydraulicNode class.
class HydraulicLink;

/**
 * @brief A Hydraulic Node in a Hydraulic Network.
 * @details A Hydraulic Node contains an unordered set of Hydraulic Links
 * connected to the Hydraulic Node, an unordered set of lateral Point Flows
 * to the Hydraulic Node, and the total energy head at the Hydraulic Node.
 *
 */
class HydraulicNode {
public:
  HydraulicNode();
  ~HydraulicNode();
  double H; /**< The energy head at this Hydraulic Node [UNITS = FT].*/
  std::unordered_set<HydraulicLink *>
      links; /**< Hydraulic Links connected to this Hydraulic Node.*/
  std::vector<double>
      point_flows; /**< The lateral point flows generating additional flow to
                         the Hydraulic Node [UNITS = CFS].*/
  double continuity(); /**< Returns the sum of the flow into the Hydraulic Node
                         minus the sum of the flow out of the Hydraulic Node
                         [UNITS = CFS].*/
};

/**
 * @brief An abstract class for a Hydraulic Link in a Hydraulic Network.
 * @details A Hydraulic Link contains the signed flow through the Hydraulic
 * Link. Concrete classes inherited from the Hydraulic Link class implement the
 * head loss calculation for that specific Hydraulic Link.
 *
 */
class HydraulicLink {
public:
  HydraulicLink();
  ~HydraulicLink();
  double Q; /**< The signed flow through the Hydraulic Link [UNITS = CFS].*/
  HydraulicNode *up_node; /**< The upstream Hydraulic Node.*/
  HydraulicNode *dn_node; /**< The downstream Hydraulic Node.*/
  std::vector<std::pair<double, double>>
      HGL; /**< The HGL profile from downstream node to upstream node.*/
  virtual double
  head_loss() = 0; /**< Calculate the head loss from the
                      downstream to the upstream Hydraulic Node [UNITS = FT].*/
  /**
   * @brief Get the maximum water surface elevation of the downstream links
   * upstream ends.
   *
   * @return double The downstream water surface elevaiton of this link.
   */
  double get_downstream_water_surface();
  /**
   * @brief Get the downstream velocity of the link.
   * @details The shape of the first downstream passage is used.
   *
   * @return double The downstream velocity [UNITS = FPS].
   */
  double get_downstream_velocity();
  void set_up_node(HydraulicNode *node);
  void set_dn_node(HydraulicNode *node);
};

/**
 * @brief A Hydraulic Component created from a network of connected Hydraulic
 * Nodes and Hydraulic Links.
 *
 */
class HydraulicComponent {
public:
  ~HydraulicComponent();
  virtual void bind(std::pair<HydraulicComponent *, unsigned int> bind_point,
                    unsigned int binding_index) = 0;
  void unbind(unsigned int binding_index);

protected:
  virtual HydraulicNode *get_binding_node(unsigned int binding_index) = 0;
  void
  bind_components(std::pair<HydraulicComponent *, unsigned int> bind_point_1,
                  std::pair<HydraulicComponent *, unsigned int> bind_point_2,
                  std::shared_ptr<HydraulicLink> link);

private:
  std::unordered_map<unsigned int,
                     std::pair<HydraulicComponent *, unsigned int>>
      binding_points;
  std::unordered_map<unsigned int, std::shared_ptr<HydraulicLink>>
      binding_links; /**< The binding links of the component.*/
};

} // namespace hazen

#endif