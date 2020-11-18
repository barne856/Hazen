#ifndef HYDRAULICNETWORK
#define HYDRAULICNETWORK

// Standard Library
#include <map>
#include <memory>
#include <mutex>
#include <thread>
#include <type_traits>
#include <utility>
#include <vector>

namespace hazen {
// Forward declare HydraulicLink for use as a pointer in HydraulicNode class.
class HydraulicLink;

/**
 * @brief A Hydraulic Node in a Hydraulic Network.
 * @details A Hydraulic Node contains a vector of Hydraulic Links
 * connected to the Hydraulic Node, a vector of lateral Point Flows
 * to the Hydraulic Node, and the total energy head at the Hydraulic Node.
 *
 */
class HydraulicNode : public std::enable_shared_from_this<HydraulicNode> {
public:
  double H{0.0}; /**< The energy head at this Hydraulic Node [UNITS = FT].*/
  std::vector<std::shared_ptr<HydraulicLink>>
      links; /**< Hydraulic Links connected to this Hydraulic Node.*/
  std::vector<double>
      point_flows; /**< The lateral point flows generating additional flow to
                         the Hydraulic Node [UNITS = CFS].*/
  std::vector<size_t> out_links(); /**< The indices of the outflowing links.*/
  std::vector<size_t> in_links();  /**< The indices of the inflowing links.*/
  virtual double continuity();     /**< Returns the sum of the flow into the
                             Hydraulic Node     minus the sum of the flow out of the
                             Hydraulic Node     [UNITS = CFS].*/
  void bind(std::shared_ptr<HydraulicNode>
                &node); /**< Bind this node to another node.*/
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
  double Q{0.0}; /**< The signed flow through the Hydraulic Link [UNITS = CFS].
               Positive flows from first to second node, negative flows from
               second to first node.*/
  std::vector<std::pair<double, double>>
      HGL; /**< The HGL profile from downstream node to upstream node.*/
  virtual double
  head_loss(HydraulicNode *node) = 0; /**< Calculate the head loss from the
 downstream to the upstream Hydraulic Node [UNITS = FT].*/
  /**
   * @brief Get the water surface elevation in this link at the given node.
   *
   * @param node The node to get the water surface elevation at.
   * @return double The water surface elevation in the link at the node.
   */
  virtual double get_water_surface_subcritical(HydraulicNode *node) = 0;
  HydraulicNode *antinode(HydraulicNode *node);
  template <size_t index> HydraulicNode *node() const {
    return std::get<index>(nodes).get();
  }
  template <size_t index> std::shared_ptr<HydraulicNode> &get_node() {
    return std::get<index>(nodes);
  }

private:
  std::pair<std::shared_ptr<HydraulicNode>, std::shared_ptr<HydraulicNode>>
      nodes; /**< The nodes connected to the link.*/
};

/**
 * @brief A Hydraulic Component created from a network of connected Hydraulic
 * Nodes and Hydraulic Links.
 *
 */
class HydraulicComponent {
public:
  virtual ~HydraulicComponent();
  std::vector<std::shared_ptr<HydraulicLink>> links;
  template <size_t index> HydraulicNode *node() const {
    return nodes[index].get();
  }
  template <size_t index> std::shared_ptr<HydraulicNode> &get_node() {
    return nodes[index];
  }

protected:
  std::vector<std::shared_ptr<HydraulicNode>> nodes;
};

/**
 * @brief A Hydraulic Network that is used to solve the continuity and head loss
 * equations for a network of Hydraulic Components.
 *
 */
class HydraulicNetwork {
public:
  /**
   * @brief Push a Hydraulic Component onto the stack. The component must
   * already be bound to something in the network.
   *
   * @param component
   */
  void push_component(std::shared_ptr<HydraulicComponent> component);
  /**
   * @brief Calculate the flow in - flow out for each node in the network
   * excluding outfall nodes.
   *
   * @param Q The vector of flows for the links in the network.
   * @return std::vector<double> The continuity equation for all nodes in the
   * network.
   */
  std::vector<double> network_continuity(std::vector<double> Q);
  double network_energy_head_sse(std::vector<double> Q);
  /**
   * @brief The solver for the Hydraulic Network.
   *
   * @return bool True if a solution was found, otherwise false.
   */
  bool solve();

private:
  void head_loss_branch(HydraulicNode *node);
  void continuity_branch(HydraulicNode *node);
  void compute_head_loss_branch(HydraulicLink *link, HydraulicNode *node);
  std::vector<std::shared_ptr<HydraulicComponent>> components;
  std::map<HydraulicNode *, std::vector<double>> node_energy_head;
  std::map<HydraulicNode *, double> node_continuity;
  std::vector<double> Q_solution;
  std::mutex node_energy_head_mutex{};
  std::mutex node_continuity_mutex{};
};

} // namespace hazen

#endif