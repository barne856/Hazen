#ifndef HYDRAULICNETWORK
#define HYDRAULICNETWORK

// Standard Library
#include <memory>
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
  double H; /**< Piezometric at this Hydraulic Node [UNITS = FT].*/
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
  double Q; /**< The signed flow through the Hydraulic Link [UNITS = CFS].*/
  HydraulicNode *up_node; /**< The upstream Hydraulic Node.*/
  HydraulicNode *dn_node; /**< The downstream Hydraulic Node.*/
  virtual double
  head_loss() = 0; /**< Calculate the head loss from the
                      downstream to the upstream Hydraulic Node [UNITS = FT].*/
};

/**
 * @brief A Hydraulic Component created from a network of connected Hydraulic
 * Nodes and Hydraulic Links.
 *
 */
class HydraulicComponent {
protected:
  std::unordered_set<std::shared_ptr<HydraulicNode>>
      nodes; /**< The Hydraulic Nodes used in the Hydraulic Component.*/
  std::unordered_set<std::shared_ptr<HydraulicLink>>
      links; /**< The Hydraulic Links used in the Hydraulic Component.*/
};

/**
 * @brief A Hydraulic Network of Hydraulic Components and Flows.
 * @details The Hydraulic Network class contains Hydraulic Components that
 * describe a Hydraulic System.
 *
 */
class HydraulicNetwork {
public:
  /**
   * @brief Add a Hydraulic Component to the Hydraulic Network.
   *
   * @param component The Hydraulic Component to add to the Hydraulic Network.
   */
  void add_component(std::shared_ptr<HydraulicComponent> component);

private:
  std::unordered_set<std::shared_ptr<HydraulicComponent>>
      components; /**< The Hydraulic Components in the Hydraulic Network.*/
};

} // namespace hazen

#endif