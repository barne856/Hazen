#include "Hazen/Hydraulics.hpp"
#include <algorithm>
#include <assert.h>
#include <exception>
#include <fstream>
#include <iostream>
#include <math.h>
#include <numeric>
#include <string>
#include <utility>
#include <vector>

namespace hazen {
Dimensionless reynolds(Density rho, Dynamic_Viscosity mu, Velocity V,
                       Length Dh) {
  return rho * V * Dh / mu;
}

Dimensionless froude(Velocity V, Length h) {
  if (h.val <= 0.0) {
    return Dimensionless(std::numeric_limits<double>::infinity());
  }
  if (h.val == std::numeric_limits<double>::infinity()) {
    return Dimensionless(0.0);
  }
  Dimensionless Fr{0.0};
  try {
    Fr = Dimensionless(std::abs(V.val) / std::sqrt(g.val * h.val));
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    Fr = Dimensionless(std::numeric_limits<double>::infinity());
  }
  return Fr;
}

Dimensionless darcy_friction_factor_pressure_driven(Dimensionless Re, Length Dh,
                                                    Length eps) {
  // Explicit approximation valid for all flow regimes from Bellos, Nalbantis,
  // and Tsakiris (2018).
  // https://ascelibrary.org/doi/10.1061/%28ASCE%29HY.1943-7900.0001540
  double a = 1.0 / (1.0 + pow(Re.val / 2712.0, 8.4));
  double b = 1.0 / (1.0 + pow(Re.val / (150.0 * Dh.val / eps.val), 1.8));
  double f1 = 64.0 / Re.val;
  double f2 = 0.75 * log(Re.val / 5.37);
  double f3 = 0.88 * log(3.41 * Dh.val / eps.val);
  return Dimensionless(pow(f1, a) * pow(f2, 2.0 * (a - 1.0) * b) *
                       pow(f3, 2.0 * (a - 1.0) * (1.0 - b)));
}

Dimensionless darcy_friction_factor_free_surface(Dimensionless Re, Length Dh,
                                                 Length eps) {
  // Explicit approximation valid for all flow regimes from Bellos, Nalbantis,
  // and Tsakiris (2018).
  // https://ascelibrary.org/doi/10.1061/%28ASCE%29HY.1943-7900.0001540
  double a = 1.0 / (1.0 + pow(Re.val / 678.0, 8.4));
  double b = 1.0 / (1.0 + pow(Re.val / (150.0 * Dh.val / eps.val), 1.8));
  double f1 = 24.0 / Re.val;
  double f2 = 0.86 * exp(lambert_W(1.35 * Re.val)) / Re.val;
  double f3 = 1.34 / pow(log(12.21 * Dh.val / eps.val), 2.0);
  return Dimensionless(pow(f1, a) * pow(f2, 2.0 * (1.0 - a) * b) *
                       pow(f3, (1.0 - a) * (1.0 - b)));
}

double lambert_W(double x) {
  double logx = log(x);
  double loglogx = log(logx);
  return logx - loglogx + loglogx / logx +
         (loglogx * loglogx - 2.0 * loglogx) / (2.0 * logx * logx);
}

Length normal_depth(HydraulicShape *shape, FrictionMethod *friction, Angle S,
                    Flow Q) {
  if (Q == Flow(0.0)) {
    return Length(std::numeric_limits<double>::infinity());
  }
  if (abs(friction->friction_slope(shape, Q, shape->get_shape_height())) >=
      abs(S)) {
    return Length(std::numeric_limits<double>::infinity());
  }
  // define function to use with bisection method
  std::function<Angle(Length)> objective = [shape, friction,
                                            Q](Length depth) -> Angle {
    return friction->friction_slope(shape, Q, depth);
  };
  double TOL = 0.00000001 * shape->get_shape_height().val;
  Length a = Dimensionless(0.5) * shape->get_shape_height();
  Length b = a;
  while (objective(a) - S < Angle(0.0)) {
    a = a / Dimensionless(2.0);
    if (a.val < TOL) {
      return Length(0.0);
    }
  }
  while (objective(b) - S > Angle(0.0)) {
    b = b * Dimensionless(2.0);
    if (b > shape->get_shape_height()) {
      return Length(std::numeric_limits<double>::infinity());
    }
  }
  // solve for depth when friction slope = bed slope
  return find_goal_bisection<Length, Angle>(S, a, b, TOL, objective);
}

Length critical_depth(HydraulicShape *shape, Flow Q) {
  Length TOL = Dimensionless(0.00000001) * shape->get_shape_height();
  if (Q.val == std::numeric_limits<double>::infinity()) {
    return Length(std::numeric_limits<double>::infinity());
  }
  if (Q.val <= TOL.val) {
    if (Q.val <= 0.0) {
      return Length(0.0);
    }
    return Length(power<2, 3>(Q).val);
  }

  // define function to use with bisection method
  std::function<Dimensionless(Length)> objective =
      [shape, Q](Length depth) -> Dimensionless {
    Area A = shape->flow_area(depth);
    if (A.val <= 0.0) {
      return Dimensionless(std::numeric_limits<double>::infinity());
    }
    Velocity V = Q / A;
    Length h = shape->hydraulic_depth(depth);
    return froude(V, h);
  };

  Dimensionless fa = (objective(TOL) - Dimensionless(1.0));
  Dimensionless fb =
      (objective(shape->get_shape_height()) - Dimensionless(1.0));
  if (fa * fb > Dimensionless(0.0)) {
    return Length(std::numeric_limits<double>::infinity());
  }
  Length a = TOL;
  Length b = shape->get_shape_height();
  // solve for depth when Fr = 1.0
  return find_goal_bisection<Length, Dimensionless>(Dimensionless(1.0), a, b,
                                                    TOL.val, objective);
}

std::vector<Vec<Length>> gen_alignment(Angle slope, Length reach, Length down_invert) {
  Vec<Length> dn_vert(3);
  Vec<Length> up_vert(3);
  dn_vert.elems(2) = down_invert.val;
  up_vert.elems(0) = reach.val;
  up_vert.elems(1) = 0.0;
  up_vert.elems(2) = slope.as_slope() * reach.val + dn_vert.elems(2);
  return {up_vert, dn_vert};
}

// std::vector<std::pair<std::string, std::vector<double>>>
// csv_util::gen_table(std::vector<std::pair<double, double>> values,
//                     std::string label1, std::string label2) {
//   std::vector<double> vec1{};
//   std::vector<double> vec2{};
//   for (auto [x, y] : values) {
//     vec1.push_back(x);
//     vec2.push_back(y);
//   }
//   std::vector<std::pair<std::string, std::vector<double>>> cols{};
//   cols.push_back({label1, vec1});
//   cols.push_back({label2, vec2});
//   return cols;
// }
//
// void csv_util::write_csv(
//     std::string filename,
//     std::vector<std::pair<std::string, std::vector<double>>> dataset) {
//   // Make a CSV file with one or more columns of integer values
//   // Each column of data is represented by the pair <column name, column
//   data>
//   //   as std::pair<std::string, std::vector<double>>
//   // The dataset is represented as a vector of these columns
//   // Note that all columns should be the same size
//
//   // Create an output filestream object
//   std::ofstream myFile(filename);
//
//   // Send column names to the stream
//   for (int j = 0; j < dataset.size(); ++j) {
//     myFile << dataset.at(j).first;
//     if (j != dataset.size() - 1)
//       myFile << ","; // No comma at end of line
//   }
//   myFile << "\n";
//
//   // Send data to the stream
//   for (int i = 0; i < dataset.at(0).second.size(); ++i) {
//     for (int j = 0; j < dataset.size(); ++j) {
//       myFile << dataset.at(j).second.at(i);
//       if (j != dataset.size() - 1)
//         myFile << ","; // No comma at end of line
//     }
//     myFile << "\n";
//   }
//
//   // Close the file
//   myFile.close();
// }

} // namespace hazen