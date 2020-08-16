#ifndef HYDRAULICNETWORK
#define HYDRAULICNETWORK

// Standard Library
#include <cmath>
#include <memory>
#include <unordered_set>

namespace hazen {
// forward declarations
class HydraulicLink;

/**
 * @brief A Hydraulic Node in a Hydraulic Network.
 * @details A Hydraulic Node contains an unordered set of Hydraulic Links
 * connected to the Hydraulic Node and the total energy head at the Hydraulic
 * Node.
 *
 */
struct HydraulicNode {
  float E{nanf("")}; /**<< Total head energy at this Hydraulic Node.*/
  std::unordered_set<HydraulicLink *>
      links; /**<< Hydraulic Links connected to this Hydraulic Node.*/
};

/**
 * @brief An abstract class for a Hydraulic Link in a Hydraulic Network.
 * @details A Hydraulic Link contatins the signed flow through the Hydraulic
 * Link and the upstream and downstream Hydraulic Nodes. Concrete classes
 * inherited from the Hydraulic Link class implement the head loss calculation
 * for that specific Hydraulic Link.
 *
 */
class HydraulicLink {
public:
  float Q{nanf("")}; /**<< The signed flow through the Hydraulic Link from the
              upstream Hydraulic Node to the downstream Hydraulic Node.*/
  HydraulicNode *up_node{nullptr}; /**<< The upstream Hydraulic Node.*/
  HydraulicNode *dn_node{nullptr}; /**<< The downstream Hydraulic Node.*/
  virtual float head_loss() = 0;   /**<< Calculate the head loss from the
                                      downstream to the upstream Hydraulic Node.*/
};

/**
 * @brief A Hydraulic Component created from a network of connected Hydraulic
 * Nodes and Hydraulic Links.
 *
 */
class HydraulicComponent {
protected:
  std::unordered_set<std::shared_ptr<HydraulicNode>> nodes;
  std::unordered_set<std::shared_ptr<HydraulicLink>> links;
};

/**
 * @brief A Point Flow that is applied to a Hydraulic Node in a Hydraulic
 * Network as a generator of additional flow to the Hydraulic Node.
 *
 */
struct PointFlow {
  float Q{nanf("")};   /**<< The additional inflow to a Hydraulic Node.*/
  HydraulicNode *node; /**<< The Hydraulic Node to apply the point flow to.*/
};

/**
 * @brief A spatially Varied Flow that is applied to a Hydraulic Link in a
 * Hydraulic Network as a generator of additional flow to the Hydraulic Link.
 * @details The spatially Varied Flow varies linearly from the upstream
 * Hydraulic Node to the downstream Hydraulic Node of the Hydraulic Link for
 * which it is applied.
 *
 */
struct VariedFlow {
  float Qup{nanf("")}; /**<< The additional inflow to a Hydraulic Link at the
                          upstream Node.*/
  float Qdn{nanf("")}; /**<< The additional inflow to a Hydraulic Node at the
                          downstream Node.*/
  HydraulicLink
      *link; /**<< The Hydraulic Link to apply the varied flow across.*/
};

class HydraulicNetwork {
public:
  void add_point_flow(std::shared_ptr<PointFlow> point_flow);
  void add_varied_flow(std::shared_ptr<VariedFlow> varied_flow);
  void add_component(std::shared_ptr<HydraulicComponent> component);

private:
  std::unordered_set<std::shared_ptr<PointFlow>> point_flows;
  std::unordered_set<std::shared_ptr<VariedFlow>> varied_flows;
  std::unordered_set<std::shared_ptr<HydraulicComponent>> components;
};

} // namespace hazen

#endif