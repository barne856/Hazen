#include "HeadLoss.hpp"
#include "HydraulicLinks.hpp"
#include "HydraulicNetwork.hpp"
#include "HydraulicShapes.hpp"
#include "HydraulicUtil.hpp"
#include <iostream>

using namespace hazen;

int main() {
  // Nodes
  std::shared_ptr<HydraulicNode> node1 = std::make_shared<HydraulicNode>();
  std::shared_ptr<HydraulicNode> node2 = std::make_shared<HydraulicNode>();
  std::shared_ptr<HydraulicNode> node3 = std::make_shared<HydraulicNode>();
  node1->H = -1.0;

  // Passage parameters
  std::shared_ptr<HydraulicShape> circle = std::make_shared<Circle>(1.0);
  std::shared_ptr<FrictionMethod> manning =
      std::make_shared<ManningsFriction>(0.013);
  vec3 v1 = {0.0, 0.0, 0.0};
  vec3 v2 = {100.0, 0.0, 0.1};
  vec3 v3 = {200.0, 0.0, 0.5};
  double ds = 0.01;

  // Passages
  std::shared_ptr<PassageLink> link1 =
      std::make_shared<PassageLink>(circle, manning, v2, v1, ds);
  std::shared_ptr<PassageLink> link2 =
      std::make_shared<PassageLink>(circle, manning, v3, v2, ds);
  link1->Q = 1.0;
  link1->dn_node = node1.get();
  link1->up_node = node2.get();
  link2->Q = 1.0;
  link2->dn_node = node2.get();
  link2->up_node = node3.get();

  // Link passages to nodes
  node1->links.insert(link1.get());
  node2->links.insert(link1.get());
  node2->links.insert(link2.get());
  node3->links.insert(link2.get());

  // run head loss calcs
  link1->head_loss();
  link2->head_loss();

  std::vector<double> X_vec{};
  std::vector<double> HGL_vec{};
  for(auto pair : link1->HGL)
  {
    X_vec.push_back(pair.first);
    HGL_vec.push_back(pair.second);
  }
  double L = link1->HGL[link1->HGL.size()-1].first;
  for(auto pair : link2->HGL)
  {
    X_vec.push_back(pair.first+L);
    HGL_vec.push_back(pair.second);
  }

  write_csv("PassageLinkTest.csv", {{"X", X_vec},
                           {"HGL", HGL_vec}});

  // double Q = 2.0;
  // double Cd = 0.62;
  // double invert = 3.0;
  // double h = -1.0;
  // double dy = 0.001;
  // std::vector<std::pair<double, double>> HGL{};
  // Rectangle orifice = hazen::Rectangle(3.0, 3.0, false);
  // double delta_H = opening_loss(&orifice, Cd, invert, Q, h, dy,
  // HGL, 1.0); std::cout << "Opening Loss = " << delta_H << std::endl;
  // std::cout << "HGL down = " << HGL[0].second << std::endl;
  // std::cout << "HGL up = " << HGL[1].second << std::endl;
  // std::cout << "Critical Depth = " << critical_depth(&orifice, Q) <<
  // std::endl;

  // Circle up_circ = hazen::Circle(1.0);
  // Circle dn_circ = hazen::Circle(0.5);
  // std::vector<std::pair<double, double>> HGL{};
  // double hL = transition_loss(&up_circ, &dn_circ, 0.0, 0.0, 1.0, 1.0,
  // 0.6, HGL); std::cout << "Transition Loss = " << hL << std::endl;
  // std::cout << "HGL down = " << HGL[0].second << std::endl; std::cout <<
  // "HGL up = " << HGL[1].second << std::endl;

  // double Q = 1.0;
  // double d = 1.0;
  // hazen::Circle circ = hazen::Circle(1.0);
  // double Dh = circ.hydraulic_diameter(d);
  // double Re = hazen::reynolds(hazen::rho_water_50, hazen::mu_water_50,
  //                             Q / circ.flow_area(d), Dh);
  // double eps = 0.0333333333333;
  // double f;
  //
  // // create data
  // Re = 1;
  // std::vector<double> Re_vec{};
  // std::vector<double> f_press{};
  // std::vector<double> f_free{};
  // std::vector<double> f_press_cw{};
  // std::vector<double> f_free_cw{};
  // for (int i = 0; i < 30; i++) {
  //   Re_vec.push_back(Re);
  //   f = hazen::darcy_friction_factor_pressure_driven(Re, Dh, eps);
  //   f_press.push_back(f);
  //   f = hazen::darcy_friction_factor_free_surface(Re, Dh, eps);
  //   f_free.push_back(f);
  //   std::function<double(double)> colebrook_white_pressure =
  //       [Re, Dh, eps](double f) -> double {
  //     return -2.0 * log10(eps / (3.7 * Dh) + 2.51 / (Re * sqrt(f))) -
  //            1.0 / sqrt(f);
  //   };
  //   std::function<double(double)> colebrook_white_free_surface =
  //       [Re, Dh, eps](double f) -> double {
  //     return -2.0 * log10(eps / (3.0 * Dh) + 2.51 / (Re * sqrt(f))) -
  //            1.0 / sqrt(f);
  //   };
  //   f = hazen::find_goal_inverse_quadratic(0.0, 0.01, 0.03, 0.05,
  //   0.0000000001,
  //                                          colebrook_white_pressure,
  //                                          100000);
  //   f_press_cw.push_back(f);
  //   f = hazen::find_goal_inverse_quadratic(0.0, 0.01, 0.03, 0.05,
  //   0.0000000001,
  //                                          colebrook_white_free_surface,
  //                                          100000);
  //   f_free_cw.push_back(f);
  //   Re *= 2.0;
  // }
  // hazen::write_csv("DarcyF.csv", {{"Re", Re_vec},
  //                                 {"f pressure", f_press},
  //                                 {"f free sureface", f_free},
  //                                 {"CW press", f_press_cw},
  //                                 {"CW free", f_press_cw}});
  return 0;
}