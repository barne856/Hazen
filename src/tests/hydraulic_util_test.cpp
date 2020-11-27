#include "HydraulicComponents.hpp"

using namespace hazen;

int main() {
  // Pipes connecting to nodes and outfalls.
  //
  //             Node2: Flow = 1.0 cfs
  //            /    \
  //  Passage3 /      \ Passage2
  //           \      /
  //            \    /
  //             Node1
  //               |
  //               | Passage1
  //               |
  //               |
  //            Outfall: H = 2.0 ft

  // Define Component Information
  auto one_foot_circle = std::make_shared<Circle>(1.0);
  auto manning_concrete = std::make_shared<ManningsFriction>(0.013);
  std::vector<vec3> alignment1 = {{0.0, 100.0, 1.0}, {0.0, 0.0, 0.0}};
  std::vector<vec3> alignment2 = {{0.0, 200.0, 2.4}, {0.0, 100.0, 1.0}};
  std::vector<vec3> alignment3 = {{0.0, 500.0, 2.4}, {0.0, 100.0, 1.0}};
  double starting_H = 2.0; // ft
  double total_Q = 1.0;    // cfs

  // Define Components
  auto outf = std::make_shared<Outfall>(starting_H);
  auto node1 = std::make_shared<Node>();
  auto node2 = std::make_shared<Node>();
  auto passage1 =
      std::make_shared<Passage>(one_foot_circle, manning_concrete, alignment1);
  auto passage2 =
      std::make_shared<Passage>(one_foot_circle, manning_concrete, alignment2);
  auto passage3 =
      std::make_shared<Passage>(one_foot_circle, manning_concrete, alignment3);

  // Bind Components
  passage1->node<Passage::NODE1>()->bind(node1->get_node<Node::NODE>());
  passage1->node<Passage::NODE2>()->bind(outf->get_node<Outfall::NODE>());
  passage2->node<Passage::NODE1>()->bind(node2->get_node<Node::NODE>());
  passage2->node<Passage::NODE2>()->bind(node1->get_node<Node::NODE>());
  passage3->node<Passage::NODE1>()->bind(node2->get_node<Node::NODE>());
  passage3->node<Passage::NODE2>()->bind(node1->get_node<Node::NODE>());

  // Add Flows
  node2->add_flow(total_Q);

  // Generate Network and Push Components
  HydraulicNetwork network{};
  network.push_component(outf);
  network.push_component(node1);
  network.push_component(node2);
  network.push_component(passage1);
  network.push_component(passage2);
  network.push_component(passage3);

  // Solve for the flow split and HGL
  network.solve();

  // return zero for success
  return 0.0;
}