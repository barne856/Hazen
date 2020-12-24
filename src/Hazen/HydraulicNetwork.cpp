#include "Hazen/HydraulicComponents.hpp"
#include "Hazen/HydraulicLinks.hpp"
#include <algorithm>
#include <chrono>
#include <iostream>
#include <iterator>
#include <numeric>
#include <set>

namespace hazen {

std::vector<size_t> HydraulicNode::out_links() const {
  std::vector<size_t> indices{};
  for (size_t i = 0; i < links.size(); i++) {
    if ((links[i]->node(0) == this && links[i]->Q > Flow(0.0)) ||
        (links[i]->node(1) == this && links[i]->Q < Flow(0.0))) {
      indices.push_back(i);
    }
  }
  return indices;
}
std::vector<size_t> HydraulicNode::in_links() const {
  std::vector<size_t> indices{};
  for (size_t i = 0; i < links.size(); i++) {
    if ((links[i]->node(0) == this && links[i]->Q <= Flow(0.0)) ||
        (links[i]->node(1) == this && links[i]->Q >= Flow(0.0))) {
      indices.push_back(i);
    }
  }
  return indices;
}
Flow HydraulicNode::continuity() {
  Flow Q{0.0};
  for (const auto &link : links) {
    if (link->node(1) == this) {
      Q += link->Q;
    }
    if (link->node(0) == this) {
      Q -= link->Q;
    }
  }
  for (const auto &flow : point_flows) {
    Q += flow;
  }
  return Q;
}

void HydraulicNode::bind(std::shared_ptr<HydraulicNode> &node) {
  if (node.get() == this) {
    // if nodes are the same, do nothing.
    return;
  } else {
    // add links to this node, and redirect them to this node.
    for (auto link : node->links) {
      if (link->node(0) == node.get()) {
        link->get_node(0) = shared_from_this();
        links.push_back(link);
      } else {
        link->get_node(1) = shared_from_this();
        links.push_back(link);
      }
    }
    // add flows to this node.
    for (auto flow : node->point_flows) {
      point_flows.push_back(flow);
    }
  }
  // set the node to this.
  node = shared_from_this();
}

HydraulicNode *HydraulicLink::antinode(HydraulicNode *node) const {
  if (this->node(0) == node) {
    return this->node(1);
  }
  if (this->node(1) == node) {
    return this->node(0);
  }
  return nullptr;
}

HydraulicComponent::~HydraulicComponent(){};

void HydraulicNetwork::push_component(
    std::shared_ptr<HydraulicComponent> component) {
  for (const auto &link : component->links) {
    link_head.try_emplace(link.get(), Length(0.0));
  }
  for (const auto &node : component->nodes) {
    node_head.try_emplace(node.get(), std::vector<Length>());
  }
  components.emplace(component);
}

void HydraulicNetwork::compute_branch_head_loss(HydraulicNode *node) {
  if (node_head[node].size() <= 1) {
    for (size_t i : node->in_links()) {
      auto link = node->links[i].get();
      auto up_node = link->antinode(node);
      Length H = link->head_loss(node);
      node_head[up_node].push_back(H);
      link_head[link] = H;
      link_upstream_node.try_emplace(link, up_node);
      compute_branch_head_loss(up_node);
    }
  }
}

void HydraulicNetwork::head_loss(const Vec<Flow> &Q) {
  // set the flow for each link
  size_t i = 0;
  for (auto &[link, _] : link_head) {
    link->Q = Q(i);
    i++;
  }

  // empty the node head map values
  for (auto &[_, head] : node_head) {
    head = {};
  }

  // clear link upstream nodes
  link_upstream_node.clear();

  // for all constant head components, initialize head and compute upstream
  // branch head losses
  for (auto &comp : components) {
    if (auto &ch_comp = std::dynamic_pointer_cast<ConstantHeadNode>(comp)) {
      auto node = ch_comp->nodes[0].get();
      node->H = ch_comp->H;
      node_head[node].push_back(ch_comp->H);
      compute_branch_head_loss(node);
    }
  }
  // if some links were not calculated, assume their upstream nodes.
  for (const auto &[link, _] : link_head) {
    link_upstream_node.try_emplace(link, link->Q >= 0.0_cfs ? link->node(0)
                                                            : link->node(1));
  }
}

Dimensionless HydraulicNetwork::Lagrangian(const Vec<Dimensionless> &x) {
  Vec<Flow> Q(link_head.size());
  for (size_t i = 0; i < Q.size(); i++) {
    Q.elems(i) = x.elems(i);
  }
  head_loss(Q);
  auto continuity_error = flow_continuity_error();
  auto head_residual = energy_head_residual();
  Vec<Dimensionless> lambda(continuity_error.size());
  size_t j = Q.size();
  for (size_t i = 0; i < continuity_error.size(); i++) {
    lambda.elems(i) = x.elems(i + j);
  }

  return Dimensionless(squared_norm(head_residual).val -
                       dot_product(lambda, continuity_error).val);
}

Vec<Dimensionless> HydraulicNetwork::objective(const Vec<Dimensionless> &x0) {

  double TOL = 1e-4;
  Vec<Dimensionless> result(x0.size());
  auto f0 = Lagrangian(x0);
  for (size_t i = 0; i < result.size(); i++) {
    auto x1 = x0;
    x1.elems(i) += TOL;
    auto f1 = Lagrangian(x1);
    result.elems(i) = (f1.val - f0.val) / TOL;
  }

  return result;

  //  Vec<Dimensionless> result(x0.size());
  //  Vec<Flow> Q(link_head.size());
  //  for (size_t i = 0; i < Q.size(); i++) {
  //    Q.elems(i) = x0.elems(i);
  //  }
  //  head_loss(Q);
  //  auto continuity_error_map = flow_continuity_error_map();
  //  auto head_residual_map = energy_head_residual_map();
  //  Vec<Dimensionless> lambda(continuity_error_map.size());
  //  size_t j = Q.size();
  //  for (size_t i = 0; i < continuity_error_map.size(); i++) {
  //    lambda.elems(i) = x0.elems(i + j);
  //  }
  //
  //  size_t i = 0;
  //  // add f - (lambda_dn - lambda_up) to result
  //  for (const auto &[link, h] : head_residual_map) {
  //    auto up_node = link_upstream_node[link];
  //    auto dn_node = link->antinode(up_node);
  //    size_t up_index = std::distance(continuity_error_map.begin(),
  //                                    continuity_error_map.find(up_node));
  //    size_t dn_index = std::distance(continuity_error_map.begin(),
  //                                    continuity_error_map.find(dn_node));
  //    double l_up =
  //        up_index < continuity_error_map.size() ? lambda.elems(up_index) :
  //        0.0;
  //    double l_dn =
  //        dn_index < continuity_error_map.size() ? lambda.elems(dn_index) :
  //        0.0;
  //    result.elems(i) = h.val - (l_dn - l_up);
  //    i++;
  //  }
  //  // add continuity to result
  //  for (const auto &[node, c] : continuity_error_map) {
  //    result.elems(i) = c.val;
  //    i++;
  //  }
  //
  //  return result;
}

Mat<Dimensionless> HydraulicNetwork::Jacobian(const Vec<Dimensionless> &x,
                                              Dimensionless dx) {
  size_t n = x.size();
  Vec<Dimensionless> zero(n);
  Mat<Dimensionless> result(n, n);
  Vec<Dimensionless> f1 = objective(x);
  for (int i = 0; i < result.n_cols(); i++) {
    Vec<Dimensionless> x2 = x;
    x2.elems(i) += dx.val;
    Vec<Dimensionless> f2 = objective(x2);
    result.elems.col(i) = ((f2 - f1) / dx).elems;
    i++;
  }
  return result;
}

Vec<Flow> HydraulicNetwork::flow_continuity_error() {
  std::vector<Flow> temp_continuity_vec{};
  for (const auto &[node, _] : node_head) {
    if (node->is_constant_head) {
      // continuity_error.elems(i) = 0.0;
      continue;
    } else {
      temp_continuity_vec.push_back(node->continuity());
    }
  }
  Vec<Flow> continuity_error(temp_continuity_vec.size());
  for (size_t i = 0; i < temp_continuity_vec.size(); i++) {
    continuity_error.elems(i) = temp_continuity_vec.at(i).val;
  }
  return continuity_error;
}

std::map<HydraulicNode *, Flow> HydraulicNetwork::flow_continuity_error_map() {
  std::map<HydraulicNode *, Flow> continuity_error_map{};
  for (const auto &[node, _] : node_head) {
    if (node->is_constant_head) {
      // continuity_error_map.try_emplace(node, 0.0_cfs);
      continue;
    } else {
      continuity_error_map.try_emplace(node, node->continuity());
    }
  }
  return continuity_error_map;
}

Vec<Length> HydraulicNetwork::energy_head_residual() {
  Vec<Length> head_residual(link_head.size());
  size_t i = 0;
  for (const auto &[link, head] : link_head) {
    auto node = link_upstream_node[link];
    auto H_vec = node_head[node];
    double H_mean = std::accumulate(H_vec.begin(), H_vec.end(), 0.0_m).val /
                    static_cast<double>(H_vec.size());
    Length res = H_vec.size() <= 1 ? Length(0.0) : Length(H_mean) - head;
    head_residual.elems(i) = res.val;
    i++;
  }
  return head_residual;
}

std::map<HydraulicLink *, Length> HydraulicNetwork::energy_head_residual_map() {
  std::map<HydraulicLink *, Length> head_residual_map{};
  size_t i = 0;
  for (const auto &[link, head] : link_head) {
    auto node = link_upstream_node[link];
    auto H_vec = node_head[node];
    double H_mean = std::accumulate(H_vec.begin(), H_vec.end(), 0.0_m).val /
                    static_cast<double>(H_vec.size());
    Length res = H_vec.size() <= 1 ? Length(0.0) : Length(H_mean) - head;
    head_residual_map.try_emplace(link, res);
    i++;
  }
  return head_residual_map;
}

void HydraulicNetwork::solve() {
  // iteration parameters
  double TOL = 1e-2;
  int MAX_ITER = 5000;
  std::cout << "Solution Begin --------------------" << std::endl;
  auto start = std::chrono::high_resolution_clock::now();

  // find least squares solution
  auto x1 = initialize_solution();
  Vec<Dimensionless> x0(x1.size());
  Vec<Dimensionless> F0 = objective(x0);
  Vec<Dimensionless> F1 = objective(x1);
  std::cout << "F =\n" << F1 << std::endl;
  x1.elems(1) += TOL;
  F1 = objective(x1);
  std::cout << "F =\n" << F1 << std::endl;
  // Broyden's method
  double convergence = 2 * TOL;
  int iters = 0;
  // iterate until convergence
  while (iters < MAX_ITER && TOL < convergence) {
    // J += ((dF - J * dx) / squared_norm(dx)) * transpose(dx);
    auto J = Jacobian(x1, Dimensionless(TOL));
    auto dx = hazen::solve_least_squares(J, -F1);
    x1 += dx;
    x0 = x1;
    F0 = F1;
    F1 = objective(x1);

    // std::cout << "J =\n" << J << std::endl;

    // affine invariant stopping criteria (works with any scale of problem)
    // convergence = dot_product(abs(F1), abs(dQ)).val;

    // non-affine invariant stopping criteria
    convergence = norm(F1).val;
    // std::cout << F1 << std::endl;
    std::cout << iters << ": " << convergence << std::endl;
    iters++;
  }
  if (iters == MAX_ITER) {
    x1 = Dimensionless(nan(""));
  }
  auto stop = std::chrono::high_resolution_clock::now();
  std::cout << "Solution End --------------------" << std::endl;
  auto duration =
      std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
  std::cout << "Solution time: " << duration.count() << " ms" << std::endl;
  std::cout << "Iterations = " << iters << std::endl;
  std::cout << "x =\n" << x1 << std::endl;
  std::cout << "Objective =\n" << F1 << std::endl;
  std::cout << "Convergence = " << convergence << std::endl;
}

} // namespace hazen