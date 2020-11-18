#include "HydraulicComponents.hpp"
#include <utility>

namespace hazen {
// Node ------------------------------------------------------------------------
Node::Node() { nodes.push_back(std::make_shared<HydraulicNode>()); }
void Node::add_flow(double Q) { nodes[0]->point_flows.push_back(Q); }

// Outfall ---------------------------------------------------------------------
Outfall::Outfall(double H) : H(H) {
  nodes.push_back(std::make_shared<HydraulicNode>());
}

// Opening ---------------------------------------------------------------------
Opening::Opening(std::shared_ptr<HydraulicShape> opening_shape, double Cd,
                 double invert) {
  auto opening_link =
      std::make_shared<OpeningLink>(std::move(opening_shape), Cd, invert);
  auto up_node = std::make_shared<HydraulicNode>();
  auto dn_node = std::make_shared<HydraulicNode>();
  opening_link->get_node<Opening::NODE1>() = up_node;
  opening_link->get_node<Opening::NODE2>() = dn_node;
  up_node->links.push_back(opening_link);
  dn_node->links.push_back(opening_link);
  nodes.push_back(dn_node);
  nodes.push_back(up_node);
  links.push_back(std::move(opening_link));
}

// Passage ---------------------------------------------------------------------
Passage::Passage(std::shared_ptr<HydraulicShape> cross_section_shape,
                 std::shared_ptr<FrictionMethod> friction_method,
                 std::vector<vec3> alignment) {
  // generate passage links
  for (int i = 0; i < alignment.size() - 1; i++) {
    auto passage_link = std::make_shared<PassageLink>(
        cross_section_shape, friction_method,
        std::make_pair(alignment[i], alignment[i + 1]));
    links.push_back(std::move(passage_link));
  }
  // generate nodes
  for (int i = 0; i < alignment.size(); i++) {
    nodes.push_back(std::make_shared<HydraulicNode>());
  }
  // set nodes of links
  for (int i = 0; i < links.size(); i++) {
    links[i]->get_node<Passage::NODE1>() = nodes[i];
    links[i]->get_node<Passage::NODE2>() = nodes[i + 1];
  }
  // set links of nodes
  nodes[0]->links.push_back(links[0]);
  for (int i = 1; i < nodes.size() - 1; i++) {
    nodes[i]->links.push_back(links[i - 1]);
    nodes[i]->links.push_back(links[i]);
  }
  nodes.back()->links.push_back(links.back());

  // make the last node the second for indexing
  std::iter_swap(nodes.begin() + nodes.size() - 1, nodes.begin() + 1);
}

} // namespace hazen
