#include "HeadLoss.hpp"
#include <fstream>
#include <iostream>
#include <string>
#include <utility> // std::pair
#include <vector>

using namespace hazen;

void write_csv(
    std::string filename,
    std::vector<std::pair<std::string, std::vector<double>>> dataset) {
  // Make a CSV file with one or more columns of integer values
  // Each column of data is represented by the pair <column name, column data>
  //   as std::pair<std::string, std::vector<double>>
  // The dataset is represented as a vector of these columns
  // Note that all columns should be the same size

  // Create an output filestream object
  std::ofstream myFile(filename);

  // Send column names to the stream
  for (int j = 0; j < dataset.size(); ++j) {
    myFile << dataset.at(j).first;
    if (j != dataset.size() - 1)
      myFile << ","; // No comma at end of line
  }
  myFile << "\n";

  // Send data to the stream
  for (int i = 0; i < dataset.at(0).second.size(); ++i) {
    for (int j = 0; j < dataset.size(); ++j) {
      myFile << dataset.at(j).second.at(i);
      if (j != dataset.size() - 1)
        myFile << ","; // No comma at end of line
    }
    myFile << "\n";
  }

  // Close the file
  myFile.close();
}

void draw_profile(HydraulicShape *shape, FrictionMethod *friction, vec3 dn_inv,
                  vec3 up_inv, double Q, double H, double dx) {
  double jump_x = nan("");
  double x = 0.0;
  double L = horizontal_length(dn_inv, up_inv);
  double S = slope(dn_inv, up_inv);
  std::vector<double> HGL{};
  std::vector<double> X{};
  std::vector<double> INV{};
  H = gvf_backwater_loss(shape, friction, up_inv, dn_inv, Q, H, dx, jump_x, 0.0,
                         0.0);
  if (isnan(jump_x)) {
    HGL.push_back(H);
    X.push_back(x);
    INV.push_back(x * S + dn_inv.z);
  }

  while (x <= L - dx) {
    H = gvf_backwater_loss(shape, friction, up_inv, dn_inv, Q, H, dx, jump_x, x,
                           x + dx);

    if (!isnan(jump_x)) {
      break;
    }
    x += dx;
    HGL.push_back(H);
    X.push_back(x);
    INV.push_back(x * S + dn_inv.z);
  }
  // check if there is a step less than dx left at the end.
  if (x < L && isnan(jump_x)) {
    // compute a partial step to reach the end of the passage.
    double dx_part = L - x;
    H = gvf_backwater_loss(shape, friction, up_inv, dn_inv, Q, H, dx_part,
                           jump_x, x, x + dx_part);
    x += dx_part;
    HGL.push_back(H);
    X.push_back(x);
    INV.push_back(x * S + dn_inv.z);
  }
  if (!isnan(jump_x)) {
    x = 0.0;
    H = gvf_frontwater_loss(shape, friction, up_inv, dn_inv, Q,
                            up_inv.z + 100.0, dx, 0.0, 0.0);
    HGL.push_back(H);
    X.push_back(L - x);
    INV.push_back(up_inv.z - x * S);
    while (x <= L - jump_x) {
      H = gvf_frontwater_loss(shape, friction, up_inv, dn_inv, Q, H, dx, x,
                              x + dx);
      x += dx;
      HGL.push_back(H);
      X.push_back(L - x);
      INV.push_back(up_inv.z - x * S);
    }
  }
  std::vector<std::pair<std::string, std::vector<double>>> vals = {
      {"X", X}, {"INV", INV}, {"HGL", HGL}};
  write_csv("HGL.csv", vals);
}

int main() {
  // Circle circ = Circle(1.0);
  // DarcyFriction darcy = DarcyFriction(0.001, rho_water_50, mu_water_50);
  // vec3 dn_inv = {0.0, 0.0, 0.0};
  // vec3 up_inv = {100.0, 0.0, 0.4};
  // double Q = 1.0;
  // double H = -0.5;
  // double dx = 1.0;
  // double S = slope(dn_inv, up_inv);
  // draw_profile(&circ, &darcy, dn_inv, up_inv, Q, H, dx);
  // std::cout << "Uniform Depth: " << normal_depth(&circ, &darcy, S, Q)
  //           << std::endl;
  // std::cout << "Critical Depth: " << critical_depth(&circ, Q) << std::endl;

  double width = 2.0;
  double height = 2.0;
  Rectangle rect = Rectangle(width, height, true);
  double Cd = 0.6228;
  double invert = 1000.0;
  double Q = 0.0;
  double H1 = 995.0;
  double dy = 0.0001;
  double percent_open = 1.0;
  double H2 = opening_loss(&rect, Cd, invert, Q, H1, dy, percent_open);
  std::cout << "Head Loss = " << H2 - H1 << std::endl;
}