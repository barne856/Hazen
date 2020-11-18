#include "HydraulicComponents.hpp"
#include "HydraulicUtil.hpp"
#include <iostream>

using namespace hazen;

int main() {
  auto outf = std::make_shared<Outfall>(2.0);
  auto node1 = std::make_shared<Node>();
  auto node2 = std::make_shared<Node>();
  auto one_foot_circle = std::make_shared<Circle>(1.0);
  auto manning_concrete = std::make_shared<ManningsFriction>(0.013);
  std::vector<vec3> alignment1 = {{0.0, 100.0, 1.0}, {0.0, 0.0, 0.0}};
  std::vector<vec3> alignment2 = {{0.0, 200.0, 2.4}, {0.0, 100.0, 1.0}};
  std::vector<vec3> alignment3 = {{0.0, 300.0, 2.4}, {0.0, 100.0, 1.0}};
  auto passage1 =
      std::make_shared<Passage>(one_foot_circle, manning_concrete, alignment1);
  auto passage2 =
      std::make_shared<Passage>(one_foot_circle, manning_concrete, alignment2);
  auto passage3 =
      std::make_shared<Passage>(one_foot_circle, manning_concrete, alignment3);
  passage1->node<Passage::NODE1>()->bind(node1->get_node<Node::NODE>());
  passage1->node<Passage::NODE2>()->bind(outf->get_node<Outfall::NODE>());
  passage2->node<Passage::NODE1>()->bind(node2->get_node<Node::NODE>());
  passage2->node<Passage::NODE2>()->bind(node1->get_node<Node::NODE>());
  passage3->node<Passage::NODE1>()->bind(node2->get_node<Node::NODE>());
  passage3->node<Passage::NODE2>()->bind(node1->get_node<Node::NODE>());
  node2->add_flow(1.0);
  HydraulicNetwork network{};
  network.push_component(outf);
  network.push_component(node1);
  network.push_component(node2);
  network.push_component(passage1);
  network.push_component(passage2);
  network.push_component(passage3);
  network.solve();
  write_csv("./output_table_1.csv",
            gen_table(passage1->links[0]->HGL, "x", "h"));
  write_csv("./output_table_3.csv",
            gen_table(passage3->links[0]->HGL, "x", "h"));
  std::cout << node1->node<0>()->H << std::endl;
  std::cout << node2->node<0>()->H << std::endl;
  return 0.0;
}