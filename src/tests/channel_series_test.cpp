#include "Hazen/HydraulicComponents.hpp"

using namespace hazen;

int main() {

  // Define Network Information
  Flow Q_total = Flow::CFS(1.0); // flow at node F
  Length H_init = 5.0_ft;        // energy head at node F

  // Define Component Information
  auto shape_1 = gen_ref<Circle>(12.0_in);
  auto friction_1 = gen_ref<ManningsFriction>(0.013_pure);
  Vec<Length> v1(3);
  Vec<Length> v2(3);
  Vec<Length> v3(3);
  Vec<Length> v4(3);
  v2.elems(0) = Length::Feet(200.0).val;
  v3.elems(0) = Length::Feet(200.0).val;
  v3.elems(2) = Length::Feet(20.0).val;
  v4.elems(0) = Length::Feet(400.0).val;
  v4.elems(2) = Length::Feet(10.0).val;
  Alignment alignment_1 = {v4, v3, v2, v1};
  std::cout << v1 << std::endl;
  std::cout << v2 << std::endl;
  std::cout << v3 << std::endl;
  std::cout << v4 << std::endl;

  // Define Components
  auto node_A = gen_ref<VariableHeadNode>();
  auto node_B = gen_ref<ConstantHeadNode>(H_init);
  auto channel_1 = gen_ref<Passage>(shape_1, friction_1, alignment_1);

  // make channel_1 closed to atmosphere
  for (auto &node : channel_1->nodes) {
    node->is_open_to_atmosphere = false;
  }

  // Bind Components
  node_A->bind(channel_1->nodes[0]);
  node_B->bind(channel_1->nodes[1]);

  // Add Flows
  node_A->add_flow(Q_total);

  // Generate Network and Push Components
  HydraulicNetwork network{};
  network.push_component(node_A);
  network.push_component(node_B);
  network.push_component(channel_1);

  // Solve for the flow split and HGL
  network.solve();

  csv_util::write_csv("./HGL_profile1.csv",
                      csv_util::gen_table(channel_1->links[0]->HGL, "X", "HGL",
                                          csv_util::CSV_UNITS::FEET));
  csv_util::write_csv("./HGL_profile2.csv",
                      csv_util::gen_table(channel_1->links[1]->HGL, "X", "HGL",
                                          csv_util::CSV_UNITS::FEET));
  csv_util::write_csv("./HGL_profile3.csv",
                      csv_util::gen_table(channel_1->links[2]->HGL, "X", "HGL",
                                          csv_util::CSV_UNITS::FEET));

  // return zero for success
  return 0.0;
}