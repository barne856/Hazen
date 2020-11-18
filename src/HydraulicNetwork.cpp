#include "HydraulicNetwork.hpp"
#include "HydraulicComponents.hpp"
#include "HydraulicLinks.hpp"
#include "HydraulicUtil.hpp"
#include <algorithm>
#include <iostream>
#include <numeric>

namespace hazen {

std::vector<size_t> HydraulicNode::out_links() {
  std::vector<size_t> indices{};
  for (size_t i = 0; i < links.size(); i++) {
    if ((links[i]->node<0>() == this && links[i]->Q > 0.0) ||
        (links[i]->node<1>() == this && links[i]->Q < 0.0)) {
      indices.push_back(i);
    }
  }
  return indices;
}
std::vector<size_t> HydraulicNode::in_links() {
  std::vector<size_t> indices{};
  for (size_t i = 0; i < links.size(); i++) {
    if ((links[i]->node<0>() == this && links[i]->Q <= 0.0) ||
        (links[i]->node<1>() == this && links[i]->Q >= 0.0)) {
      indices.push_back(i);
    }
  }
  return indices;
}
double HydraulicNode::continuity() {
  double Q = 0.0;
  for (const auto &link : links) {
    if (link->node<1>() == this) {
      Q += link->Q;
    }
    if (link->node<0>() == this) {
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
      if (link->node<0>() == node.get()) {
        link->get_node<0>() = shared_from_this();
        links.push_back(link);
      } else {
        link->get_node<1>() = shared_from_this();
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

HydraulicNode *HydraulicLink::antinode(HydraulicNode *node) {
  if (this->node<0>() == node) {
    return this->node<1>();
  }
  if (this->node<1>() == node) {
    return this->node<0>();
  }
  return nullptr;
}

HydraulicComponent::~HydraulicComponent(){};

void HydraulicNetwork::push_component(
    std::shared_ptr<HydraulicComponent> component) {
  for (int i = 0; i < component->links.size(); i++) {
    Q_solution.push_back(0.0);
  }
  components.push_back(component);
}

double HydraulicNetwork::network_energy_head_sse(std::vector<double> Q) {
  // clear energy head maps.
  node_energy_head.clear();
  // get all outfalls and add them to the node energy head map
  std::vector<HydraulicNode *> outfalls{};
  for (auto &comp : components) {
    if (auto &outfall = std::dynamic_pointer_cast<Outfall>(comp)) {
      outfall->node<Outfall::NODE>()->H = outfall->H;
      outfalls.push_back(outfall->node<Outfall::NODE>());
      node_energy_head.try_emplace(outfall->node<Outfall::NODE>(), outfall->H);
    }
  }
  // set the flow of every link
  int k = 0;
  for (auto &comp : components) {
    for (auto &link : comp->links) {
      link->Q = Q[k];
      k++;
    }
  }
  // compute energy head at all nodes
  std::vector<std::thread> outfall_threads{};
  for (auto outfall : outfalls) {
    std::thread t{[=] { this->head_loss_branch(outfall); }};
    outfall_threads.push_back(std::move(t));
  }
  for (auto &t : outfall_threads) {
    t.join();
  }
  // Return sum of squared errors from the mean energy head for each node
  double sse = 0.0;
  for (const auto &entry : node_energy_head) {
    if (std::any_of(
            outfalls.begin(), outfalls.end(),
            [=](HydraulicNode *node) -> bool { return node == entry.first; })) {
      continue;
    }
    const auto &H_vec = entry.second;
    double H_mean = std::accumulate(H_vec.begin(), H_vec.end(), 0.0) /
                    static_cast<double>(H_vec.size());
    for (const auto &H : H_vec) {
      sse += (H_mean - H) * (H_mean - H);
    }
  }
  return sse;
}

void HydraulicNetwork::head_loss_branch(HydraulicNode *node) {
  auto out_links = node->out_links();
  auto in_links = node->in_links();
  node_energy_head_mutex.lock();
  if (node_energy_head.find(node) == node_energy_head.end() ||
      node_energy_head[node].size() < out_links.size()) {
    node_energy_head_mutex.unlock();
    return;
  }
  node_energy_head_mutex.unlock();
  std::vector<std::thread> branch_threads{};
  for (auto link : in_links) {
    std::thread t{
        [=] { this->compute_head_loss_branch(node->links[link].get(), node); }};
    branch_threads.push_back(std::move(t));
  }
  for (auto &t : branch_threads) {
    t.join();
  }
}

void HydraulicNetwork::compute_head_loss_branch(HydraulicLink *link,
                                                HydraulicNode *node) {
  double H = link->head_loss(node);
  HydraulicNode *up_node = link->antinode(node);
  node_energy_head_mutex.lock();
  up_node->H = H;
  std::vector<double> temp{};
  node_energy_head.try_emplace(up_node, temp);
  node_energy_head[up_node].push_back(H);
  node_energy_head_mutex.unlock();
  head_loss_branch(up_node);
}

std::vector<double>
HydraulicNetwork::network_continuity(std::vector<double> Q) {
  // clear node continuity map
  node_continuity.clear();
  // Assign flow to all links
  int k = 0;
  for (auto &comp : components) {
    for (auto &link : comp->links) {
      link->Q = Q[k];
      k++;
    }
  }
  // get all outfalls
  std::vector<HydraulicNode *> outfalls{};
  for (auto &comp : components) {
    if (auto &outfall = std::dynamic_pointer_cast<Outfall>(comp)) {
      outfalls.push_back(outfall->node<Outfall::NODE>());
    }
  }
  // compute continuity in parallel
  std::vector<std::thread> outfall_threads{};
  for (auto outfall : outfalls) {
    std::thread t{[=] { this->continuity_branch(outfall); }};
    outfall_threads.push_back(std::move(t));
  }
  for (auto &t : outfall_threads) {
    t.join();
  }
  // report the results
  std::vector<double> result{};
  for (const auto &entry : node_continuity) {
    if (std::any_of(
            outfalls.begin(), outfalls.end(),
            [=](HydraulicNode *node) -> bool { return node == entry.first; })) {
      continue;
    }
    result.push_back(entry.second);
  }
  return result;
}

void HydraulicNetwork::continuity_branch(HydraulicNode *node) {
  node_continuity_mutex.lock();
  if (node_continuity.find(node) == node_continuity.end()) {
    node_continuity.try_emplace(node, 0.0);
  } else {
    node_continuity_mutex.unlock();
    return;
  }
  node_continuity_mutex.unlock();
  double c = node->continuity();
  node_continuity_mutex.lock();
  node_continuity[node] = c;
  node_continuity_mutex.unlock();
  std::vector<std::thread> branch_threads{};
  for (auto link : node->links) {
    std::thread t{[=] { this->continuity_branch(link->antinode(node)); }};
    branch_threads.push_back(std::move(t));
  }
  for (auto &t : branch_threads) {
    t.join();
  }
}

bool HydraulicNetwork::solve() {
  // initial guess is the least squares solution to the network.
  std::vector<double> Q_guess = solve_least_squares_linear_system(
      Q_solution, 0.000000001,
      std::bind(&HydraulicNetwork::network_continuity, this,
                std::placeholders::_1));
  // use method of lagrange multipliers to optimize solution.
  Q_solution = solve_constrained_problem_lagrange(
      Q_guess,
      std::bind(&HydraulicNetwork::network_energy_head_sse, this,
                std::placeholders::_1),
      std::bind(&HydraulicNetwork::network_continuity, this,
                std::placeholders::_1),
      0.000001, 1000);
  std::cout << network_energy_head_sse(Q_solution) << std::endl;
  return true;
}

} // namespace hazen