#include "Hazen/HydraulicComponents.hpp"

#include <iostream>
#include <stdlib.h> /* srand, rand */
#include <time.h>   /* time */

using namespace hazen;

int main() {
  // srand(time(NULL));
  // int N = 20;
  // Mat<Dimensionless> A(N, N);
  // for (int i = 0; i < N; i++) {
  //   for (int j = 0; j < N; j++) {
  //     A(i, j) = Dimensionless(static_cast<double>(rand() % 100));
  //   }
  // }
  // std::cout << "inversion begin: " << std::endl;
  // auto B = inv(A);
  // std::cout << "inversion end: " << std::endl;

  // Define Network Information
  Flow Q_total = Flow::MGD(3.5); // flow at node F
  Length H_init = 5.1_m;         // energy head at node F

  // Define Component Information
  auto shape_1 = gen_ref<Circle>(Length::Inches(18.0));
  auto shape_2 = gen_ref<Circle>(Length::Inches(24.0));
  auto shape_3 = gen_ref<Circle>(Length::Inches(18.0));
  auto shape_4 = gen_ref<Circle>(Length::Inches(12.0));
  auto shape_5 = gen_ref<Circle>(Length::Inches(12.0));
  auto shape_6 = gen_ref<Circle>(Length::Inches(12.0));
  auto shape_7 = gen_ref<Circle>(Length::Inches(18.0));
  auto shape_8 = gen_ref<Circle>(Length::Inches(24.0));
  auto friction_1 = gen_ref<ManningsFriction>(0.013_pure);
  auto friction_2 = gen_ref<ManningsFriction>(0.013_pure);
  auto friction_3 = gen_ref<ManningsFriction>(0.012_pure);
  auto friction_4 = gen_ref<ManningsFriction>(0.014_pure);
  auto friction_5 = gen_ref<ManningsFriction>(0.013_pure);
  auto friction_6 = gen_ref<ManningsFriction>(0.013_pure);
  auto friction_7 = gen_ref<ManningsFriction>(0.014_pure);
  auto friction_8 = gen_ref<ManningsFriction>(0.014_pure);
  auto alignment_1 = gen_alignment(Angle::Slope(0.0005), 200.0_m);
  auto alignment_2 = gen_alignment(Angle::Slope(0.0005), 200.0_m);
  auto alignment_3 = gen_alignment(Angle::Slope(0.0005), 200.0_m);
  auto alignment_4 = gen_alignment(Angle::Slope(0.0005), 100.0_m);
  auto alignment_5 = gen_alignment(Angle::Slope(0.0005), 100.0_m);
  auto alignment_6 = gen_alignment(Angle::Slope(0.0005), 100.0_m);
  auto alignment_7 = gen_alignment(Angle::Slope(0.0005), 100.0_m);
  auto alignment_8 = gen_alignment(Angle::Slope(0.0005), 300.0_m);

  // Define Components
  auto node_F = gen_ref<ConstantHeadNode>(H_init);
  auto node_A = gen_ref<VariableHeadNode>();
  auto node_B = gen_ref<VariableHeadNode>();
  auto node_C = gen_ref<VariableHeadNode>();
  auto node_D = gen_ref<VariableHeadNode>();
  auto node_E = gen_ref<VariableHeadNode>();
  auto channel_1 = gen_ref<Passage>(shape_1, friction_1, alignment_1);
  auto channel_2 = gen_ref<Passage>(shape_2, friction_2, alignment_2);
  auto channel_3 = gen_ref<Passage>(shape_3, friction_3, alignment_3);
  auto channel_4 = gen_ref<Passage>(shape_4, friction_4, alignment_4);
  auto channel_5 = gen_ref<Passage>(shape_5, friction_5, alignment_5);
  auto channel_6 = gen_ref<Passage>(shape_6, friction_6, alignment_6);
  auto channel_7 = gen_ref<Passage>(shape_7, friction_7, alignment_7);
  auto channel_8 = gen_ref<Passage>(shape_8, friction_8, alignment_8);

  // Bind Components
  node_A->bind(channel_1->nodes[0]);
  node_B->bind(channel_1->nodes[1]);
  node_A->bind(channel_2->nodes[0]);
  node_D->bind(channel_2->nodes[1]);
  node_B->bind(channel_3->nodes[0]);
  node_E->bind(channel_3->nodes[1]);
  node_B->bind(channel_4->nodes[0]);
  node_C->bind(channel_4->nodes[1]);
  node_C->bind(channel_5->nodes[0]);
  node_E->bind(channel_5->nodes[1]);
  node_D->bind(channel_6->nodes[0]);
  node_C->bind(channel_6->nodes[1]);
  node_E->bind(channel_7->nodes[0]);
  node_F->bind(channel_7->nodes[1]);
  node_D->bind(channel_8->nodes[0]);
  node_F->bind(channel_8->nodes[1]);

  // Add Flows
  node_A->add_flow(Q_total);

  // Generate Network and Push Components
  HydraulicNetwork network{};
  network.push_component(node_A);
  network.push_component(node_B);
  network.push_component(node_C);
  network.push_component(node_D);
  network.push_component(node_E);
  network.push_component(node_F);
  network.push_component(channel_1);
  network.push_component(channel_2);
  network.push_component(channel_3);
  network.push_component(channel_4);
  network.push_component(channel_5);
  network.push_component(channel_6);
  network.push_component(channel_7);
  network.push_component(channel_8);

  // Solve for the flow split and HGL
  network.solve();

  // return zero for success
  return 0.0;
}