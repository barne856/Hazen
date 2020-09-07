#include "HeadLoss.hpp"
#include "HydraulicLinks.hpp"
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

std::vector<std::pair<std::string, std::vector<double>>>
draw_profile(HydraulicShape *shape, FrictionMethod *friction, vec3 dn_inv,
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
  return vals;
}

bool is_steep(HydraulicShape *shape, FrictionMethod *friction, double Q,
              vec3 dn_inv, vec3 up_inv) {
  double S = slope(dn_inv, up_inv);
  double yc = critical_depth(shape, Q);
  double yn = normal_depth(shape, friction, S, Q);
  if (yn < yc) {
    return true;
  }
  return false;
}

std::vector<std::pair<std::string, std::vector<double>>>
draw_passage_profile(PassageLink *passage) {
  double jump_x = nan("");
  int num_sections = passage->invert_alignment.size() - 1;
  auto shape = passage->cross_section_shape.get();
  auto friction = passage->friction_method.get();
  double H = passage->dn_node->H;
  double dx = 1.0;
  double Q = passage->Q;
  double total_x = 0.0;
  std::vector<double> HGL{};
  std::vector<double> X{};
  std::vector<double> INV{};
  for (int i = 0; i < num_sections; i++) {
    vec3 up_inv = passage->invert_alignment[i + 1];
    vec3 dn_inv = passage->invert_alignment[i];
    double L = horizontal_length(dn_inv, up_inv);
    double S = slope(dn_inv, up_inv);
    double x = 0.0;
    if (i == 0) {
      // first point
      H = gvf_backwater_loss(shape, friction, up_inv, dn_inv, Q, H, dx, jump_x,
                             x, x);
      if (isnan(jump_x)) {
        HGL.push_back(H);
        X.push_back(x);
        INV.push_back(x * S + dn_inv.z);
      }
    }
    while (x < L) {
      if (x + dx > L) {
        H = gvf_backwater_loss(shape, friction, up_inv, dn_inv, Q, H, dx,
                               jump_x, x);
        total_x += (L - x);
        x = L;
      } else {
        H = gvf_backwater_loss(shape, friction, up_inv, dn_inv, Q, H, dx,
                               jump_x, x, x + dx);
        x += dx;
        total_x += dx;
      }
      HGL.push_back(H);
      X.push_back(total_x);
      INV.push_back(x * S + dn_inv.z);
      if (!isnan(jump_x)) {
        // jump occured, frontwater calcs from nearest mild slope upstream
        total_x += (L - x);
        int old_i = i;
        i++;
        while (i < num_sections) {
          up_inv = passage->invert_alignment[i + 1];
          dn_inv = passage->invert_alignment[i];
          S = slope(dn_inv, up_inv);
          if (is_steep(shape, friction, Q, dn_inv, up_inv)) {
            total_x += horizontal_length(dn_inv, up_inv);
            i++;
            continue;
          }
          break;
        }
        i--;
        up_inv = passage->invert_alignment[i + 1];
        dn_inv = passage->invert_alignment[i];
        S = slope(dn_inv, up_inv);
        H = critical_depth(shape, Q) * sqrt(S * S + 1) + up_inv.z;
        HGL.push_back(H);
        X.push_back(total_x);
        INV.push_back(up_inv.z);
        double x_end_of_steep_slopes = total_x;
        for (int j = i; j >= old_i; j--) {
          up_inv = passage->invert_alignment[j + 1];
          dn_inv = passage->invert_alignment[j];
          S = slope(dn_inv, up_inv);
          double section_L = horizontal_length(dn_inv, up_inv);
          x = 0.0;
          double x_final = section_L;
          if (j == old_i) {
            x_final = section_L - jump_x;
          }
          while (x < x_final) {
            double step = dx;
            if (x + dx > x_final) {
              step = x_final - x;
            }
            H = gvf_frontwater_loss(shape, friction, up_inv, dn_inv, Q, H, dx,
                                    x, x + dx);
            x += step;
            total_x -= step;
            HGL.push_back(H);
            X.push_back(total_x);
            INV.push_back(up_inv.z - x * S);
          }
        }
        total_x = x_end_of_steep_slopes;
        break;
      }
    }
  }
  std::vector<std::pair<std::string, std::vector<double>>> vals = {
      {"X", X}, {"INV", INV}, {"HGL", HGL}};
  return vals;
}

int main() {
  PassageLink pipe = PassageLink();
  std::shared_ptr<HydraulicNode> dn_node = std::make_shared<HydraulicNode>();
  std::shared_ptr<HydraulicNode> up_node = std::make_shared<HydraulicNode>();
  std::shared_ptr<Circle> shape = std::make_shared<Circle>(5.0);
  std::shared_ptr<ManningsFriction> friction =
      std::make_shared<ManningsFriction>(0.013);
  double Q = 400.0;
  dn_node->H = 0.0;
  alignment align = alignment({{0.0, 0.0, 0.0},
                               {100.0, 0.0, 0.0},
                               {200.0, 0.0, 11.0},
                               {300.0, 0.0, 11.5},
                               {400.0, 0.0, 11.5},
                               {500.0, 0.0, 16.0},
                               {600.0, 0.0, 32.0},
                               {700.0, 0.0, 32.0},
                               {800.0, 0.0, 32.0},
                               {900.0, 0.0, 32.0},
                               {1000.0, 0.0, 32.0}});
  pipe.set_cross_section(std::move(shape));
  pipe.set_friction_method(std::move(friction));
  pipe.invert_alignment = align;
  pipe.Q = Q;
  pipe.up_node = up_node.get();
  pipe.dn_node = dn_node.get();
  auto vals = draw_passage_profile(&pipe);
  write_csv("HGL.csv", vals);

  // double width = 2.0;
  // double height = 2.0;
  // Rectangle rect = Rectangle(width, height, true);
  // double Cd = 0.6228;
  // double invert = 1000.0;
  // double Q = 0.0;
  // double H1 = 995.0;
  // double dy = 0.0001;
  // double percent_open = 1.0;
  // double H2 = opening_loss(&rect, Cd, invert, Q, H1, dy, percent_open);
  // std::cout << "Head Loss = " << H2 - H1 << std::endl;
}