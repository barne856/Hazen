#ifndef HAZEN_HYDRAULIC_NETWORK
#define HAZEN_HYDRAULIC_NETWORK

#include "Hazen/Hydraulics.hpp"

// Standard Library
#include <map>
#include <memory>
#include <mutex>
#include <thread>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>
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
  Length H{0.0}; /**< The energy head at this Hydraulic Node [UNITS = FT].*/
  std::vector<std::shared_ptr<HydraulicLink>>
      links; /**< Hydraulic Links connected to this Hydraulic Node.*/
  std::vector<Flow>
      point_flows; /**< The lateral point flows generating additional flow to
                         the Hydraulic Node [UNITS = CFS].*/
  std::vector<size_t>
  out_links() const; /**< The indices of the outflowing links.*/
  std::vector<size_t>
  in_links() const;          /**< The indices of the inflowing links.*/
  virtual Flow continuity(); /**< Returns the sum of the flow into the
                         Hydraulic Node     minus the sum of the flow out of the
                         Hydraulic Node     [UNITS = CFS].*/
  void bind(std::shared_ptr<HydraulicNode>
                &node); /**< Bind this node to another node.*/
  bool is_constant_head = false;
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
  Flow Q{0.0}; /**< The signed flow through the Hydraulic Link [UNITS = CFS].
               Positive flows from first to second node, negative flows from
               second to first node.*/
  std::vector<std::pair<Length, Length>>
      HGL; /**< The HGL profile from downstream node to upstream node.*/
  virtual Length
  head_loss(HydraulicNode *node) = 0; /**< Calculate the head loss from the
 downstream to the upstream Hydraulic Node [UNITS = FT].*/
  /**
   * @brief Get the water surface elevation in this link at the given node.
   *
   * @param node The node to get the water surface elevation at.
   * @return double The water surface elevation in the link at the node.
   */
  virtual Length get_water_surface_subcritical(HydraulicNode *node) const = 0;
  HydraulicNode *antinode(HydraulicNode *node) const;
  HydraulicNode *node(size_t index) const {
    switch (index) {
    case 0:
      return nodes.first.get();
      break;
    case 1:
      return nodes.second.get();
      break;
    default:
      return nodes.second.get();
      break;
    }
  }
  std::shared_ptr<HydraulicNode> &get_node(size_t index) {
    switch (index) {
    case 0:
      return nodes.first;
      break;
    case 1:
      return nodes.second;
      break;
    default:
      return nodes.second;
      break;
    }
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
   * @brief Push a Hydraulic Component onto the stack.
   *
   * @param component The Component to puch onto the stack,
   */
  void push_component(std::shared_ptr<HydraulicComponent> component);

  /**
   * @brief Solve the network flows split.
   *
   */
  void solve();

private:
  std::unordered_set<std::shared_ptr<HydraulicComponent>> components;
  std::map<HydraulicLink *, Length> link_head;
  std::map<HydraulicNode *, std::vector<Length>> node_head;
  std::unordered_map<HydraulicLink *, HydraulicNode *> link_upstream_node;

  /**
   * @brief Compute the head loss in a branch of the network starting from a
   * downstream node and stoping at the first previously visited upstream node
   *
   * @param node The starting node.
   */
  void compute_branch_head_loss(HydraulicNode *node);

  /**
   * @brief Compute the head loss in the network. Updates link head and node
   * head maps
   *
   * @param Q The flow to use when updating the network head loss.
   */
  void head_loss(const Vec<Flow> &Q);

  /**
   * @brief Objective function used to solve flow splits in the network.
   *
   * @param Q The flow for each link in the network.
   * @return Vec<Dimensionless> The value of the objective function. This is the
   * zero vector when a solution is reached.
   */
  Vec<Dimensionless> objective(const Vec<Dimensionless> &x0);

  /**
   * @brief Compute an approximate Jacobian Matrix of the objective function.
   *
   * @param Q1 Flow at which to approximate the Jacobian
   * @param dQ The finite difference step used to compute the Jacobian
   * @return Mat<quotient_type<Dimensionless, Flow>> The approximate Jacobian
   * Matrix for the objective function.
   */
  Mat<Dimensionless> Jacobian(const Vec<Dimensionless> &x, Dimensionless dx);

  /**
   * @brief Compute Flow Continuity Error at every node.
   *
   */
  Vec<Flow> flow_continuity_error();

  /**
   * @brief Compute Flow Continuity Error at every node.
   *
   */
  std::map<HydraulicNode *, Flow> flow_continuity_error_map();

  /**
   * @brief Compute Energy Head Residuals upstream of every link.
   *
   */
  Vec<Length> energy_head_residual();

  /**
   * @brief Compute Energy Head Residuals upstream of every link.
   *
   */
  std::map<HydraulicLink *, Length> energy_head_residual_map();

  Vec<Dimensionless> initialize_solution() {
    Vec<Flow> Q(link_head.size());
    head_loss(Q);
    auto CE = flow_continuity_error();
    Mat<Dimensionless> J(CE.size(), Q.size());
    for (size_t i = 0; i < Q.size(); i++) {
      auto Q1 = Q;
      Q1.elems(i) += 1.0;
      auto CE1 = flow_continuity_error();
      J.elems.col(i) = (CE1 - CE).elems;
    }
    // Vec<Flow> zero(CE.size());
    Q = solve_least_squares(J, -CE);

    size_t i = link_head.size();
    for (const auto &[node, _] : node_head) {
      if (!node->is_constant_head) {
        i++;
      }
    }
    Vec<Dimensionless> result(i);
    for (size_t j = 0; j < link_head.size(); j++) {
      result.elems(j) = Q.elems(j);
      if(std::abs(Q.elems(j)) < 1e-3)
      {
        result.elems(j) = 0.0;
      }
    }
    for (size_t j = link_head.size(); j < result.size(); j++) {
      result.elems(j) = 1.0;
    }
    return result;
  }
};

} // namespace hazen

#endif