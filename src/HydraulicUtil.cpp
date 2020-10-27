#include "HydraulicUtil.hpp"
#include <algorithm>
#include <assert.h>
#include <exception>
#include <fstream>
#include <iostream>
#include <math.h>
#include <numeric>
#include <string>
#include <utility> // std::pair
#include <vector>

namespace hazen {
double reynolds(double rho, double mu, double V, double Dh) {
  return rho * V * Dh / mu;
}

double froude(double V, double h) {
  if (h <= 0.0) {
    return std::numeric_limits<double>::infinity();
  }
  if (h == std::numeric_limits<double>::infinity()) {
    return 0.0;
  }
  double Fr;
  try {
    Fr = abs(V) / sqrt(g * h);
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    Fr = std::numeric_limits<double>::infinity();
  }
  return Fr;
}

double darcy_friction_factor_pressure_driven(double Re, double Dh, double eps) {
  // Explicit approximation valid for all flow regimes from Bellos, Nalbantis,
  // and Tsakiris (2018).
  // https://ascelibrary.org/doi/10.1061/%28ASCE%29HY.1943-7900.0001540
  double a = 1.0 / (1.0 + pow(Re / 2712.0, 8.4));
  double b = 1.0 / (1.0 + pow(Re / (150.0 * Dh / eps), 1.8));
  double f1 = 64.0 / Re;
  double f2 = 0.75 * log(Re / 5.37);
  double f3 = 0.88 * log(3.41 * Dh / eps);
  return pow(f1, a) * pow(f2, 2.0 * (a - 1.0) * b) *
         pow(f3, 2.0 * (a - 1.0) * (1.0 - b));
}

double darcy_friction_factor_free_surface(double Re, double Dh, double eps) {
  // Explicit approximation valid for all flow regimes from Bellos, Nalbantis,
  // and Tsakiris (2018).
  // https://ascelibrary.org/doi/10.1061/%28ASCE%29HY.1943-7900.0001540
  double a = 1.0 / (1.0 + pow(Re / 678, 8.4));
  double b = 1.0 / (1.0 + pow(Re / (150.0 * Dh / eps), 1.8));
  double f1 = 24.0 / Re;
  double f2 = 0.86 * exp(lambert_W(1.35 * Re)) / Re;
  double f3 = 1.34 / pow(log(12.21 * Dh / eps), 2.0);
  return pow(f1, a) * pow(f2, 2.0 * (1.0 - a) * b) *
         pow(f3, (1.0 - a) * (1.0 - b));
}

double lambert_W(double x) {
  double logx = log(x);
  double loglogx = log(logx);
  return logx - loglogx + loglogx / logx +
         (loglogx * loglogx - 2.0 * loglogx) / (2.0 * logx * logx);
}

double normal_depth(HydraulicShape *shape, FrictionMethod *friction, double S,
                    double Q) {
  if (Q == 0.0) {
    return std::numeric_limits<double>::infinity();
  }
  if (abs(friction->friction_slope(shape, Q, shape->get_shape_height())) >=
      abs(S)) {
    return std::numeric_limits<double>::infinity();
  }
  // define function to use with bisection method
  std::function<double(double)> objective = [shape, friction,
                                             Q](double depth) -> double {
    return friction->friction_slope(shape, Q, depth);
  };
  double TOL = 0.00000001 * shape->get_shape_height();
  double a = 0.5 * shape->get_shape_height();
  double b = a;
  while (objective(a) - S < 0.0) {
    a /= 2.0;
    if (a < TOL) {
      return 0.0;
    }
  }
  while (objective(b) - S > 0.0) {
    b *= 2;
    if (b > shape->get_shape_height()) {
      return std::numeric_limits<double>::infinity();
    }
  }
  // solve for depth when friction slope = bed slope
  return find_goal_bisection(S, a, b, TOL, objective);
}

double critical_depth(HydraulicShape *shape, double Q) {
  if (Q == std::numeric_limits<double>::infinity()) {
    return std::numeric_limits<double>::infinity();
  }
  if (Q == 0.0) {
    return 0.0;
  }
  // define function to use with bisection method
  std::function<double(double)> objective = [shape, Q](double depth) -> double {
    double A = shape->flow_area(depth);
    if (A <= 0.0) {
      return std::numeric_limits<double>::infinity();
    }
    double V = Q / A;
    double h = shape->hydraulic_depth(depth);
    return froude(V, h);
  };
  double TOL = 0.00000001 * shape->get_shape_height();
  double fa = (objective(TOL) - 1.0);
  double fb = (objective(shape->get_shape_height()) - 1.0);
  if (fa * fb > 0.0) {
    return std::numeric_limits<double>::infinity();
  }
  double a = TOL;
  double b = shape->get_shape_height();
  // solve for depth when Fr = 1.0
  return find_goal_bisection(1.0, a, b, TOL, objective);
}

double brink_depth(HydraulicShape *shape, FrictionMethod *friction, double S,
                   double Q) {
  double Dn = normal_depth(shape, friction, S, Q);
  double Dc = critical_depth(shape, Q);
  if (isnan(Dn) || isinf(Dn)) {
    return Dc;
  }
  return std::min(Dn, Dc);
}

double length(vec3 down, vec3 up) {
  return sqrt(pow(down.x - up.x, 2.0) + pow(down.y - up.y, 2.0) +
              pow(down.z - up.z, 2.0));
}

double horizontal_length(vec3 down, vec3 up) {
  return sqrt(pow(down.x - up.x, 2.0) + pow(down.y - up.y, 2.0));
}

double slope(vec3 down, vec3 up) {
  double L = horizontal_length(down, up);
  if (L == 0.0) {
    return (up.z - down.z) * std::numeric_limits<double>::infinity();
  }
  return (up.z - down.z) / (horizontal_length(down, up));
}

class no_root : public std::exception {
  virtual const char *what() const throw() {
    return "No root found before maximum iterations reached by root finding "
           "algorithm.";
  }
} nr_except;

double find_goal_secant(double goal, double x0, double x1, double TOL,
                        std::function<double(double)> objective,
                        const int MAX_ITER) {
  double result;
  try {
    double convergence = 2 * TOL;
    int iters = 0;
    double x2;
    double f_x0 = objective(x0) - goal;
    double f_x1 = objective(x1) - goal;
    while (iters < MAX_ITER && TOL < convergence) {
      x2 = x1 - f_x1 * (x1 - x0) / (f_x1 - f_x0);
      x0 = x1;
      x1 = x2;
      f_x0 = f_x1;
      f_x1 = objective(x1) - goal;
      convergence = abs(f_x1);
      iters++;
    }
    if (iters == MAX_ITER) {
      throw nr_except;
    }
    result = x2;
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    result = nan("");
  }
  return result;
}

double find_goal_bisection(double goal, double a, double b, double TOL,
                           std::function<double(double)> objective,
                           const int MAX_ITER) {
  double result;
  try {
    int i = 0;
    double c;
    double fc;
    while (i < MAX_ITER) {
      c = (a + b) / 2.0;
      fc = objective(c) - goal;
      if (fc == 0.0 || (b - c) / 2.0 < TOL) {
        // solution found
        break;
      }
      fc *(objective(a) - goal) > 0.0 ? a = c : b = c;
      i++;
    }
    if (i == MAX_ITER) {
      throw nr_except;
    }
    result = c;
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    result = nan("");
  }
  return result;
}

double find_goal_inverse_quadratic(double goal, double x0, double x1, double x2,
                                   double TOL,
                                   std::function<double(double)> objective,
                                   const int MAX_ITER) {
  double result;
  try {
    double convergence = 2 * TOL;
    int iters = 0;
    double x3;
    double f_x0 = objective(x0) - goal;
    double f_x1 = objective(x1) - goal;
    double f_x2 = objective(x2) - goal;
    while (iters < MAX_ITER && TOL < convergence) {
      double a = f_x1 * f_x2 / ((f_x0 - f_x1) * (f_x0 - f_x2));
      double b = f_x0 * f_x2 / ((f_x1 - f_x0) * (f_x1 - f_x2));
      double c = f_x0 * f_x1 / ((f_x2 - f_x0) * (f_x2 - f_x1));
      x3 = a * x0 + b * x1 + c * x2;
      x0 = x1;
      x1 = x2;
      x2 = x3;
      f_x0 = f_x1;
      f_x1 = f_x2;
      f_x2 = objective(x2) - goal;
      convergence = abs(f_x2);
      iters++;
    }
    if (iters == MAX_ITER) {
      throw nr_except;
    }
    result = x3;
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    result = nan("");
  }
  return result;
}

std::vector<double> unconstrained_solver_secant(
    std::vector<double> x0, double TOL,
    std::function<double(std::vector<double>)> objective, const int MAX_ITER) {
  // allocate space for result and fill with zeros
  size_t n = x0.size();
  std::vector<double> result{};
  result.resize(n);
  // allocate space for x1 and initialize to x0 + TOL
  std::vector<double> x1 = x0;
  for (int i = 0; i < n; i++) {
    x1[i] += TOL;
  }
  // allocate space for x2 and fill with zeros
  std::vector<double> x2;
  x2.resize(n);

  // iterate until convergence
  try {
    double convergence = 2 * TOL;
    int iters = 0;
    double f_x0 = objective(x0);
    double f_x1 = objective(x1);
    while (iters < MAX_ITER && TOL < convergence) {
      double t = f_x1 / (f_x1 - f_x0);
      for (int i = 0; i < n; i++) {
        x2[i] = x1[i] - t * (x1[i] - x0[i]);
      }
      x0 = x1;
      x1 = x2;
      f_x0 = f_x1;
      f_x1 = objective(x1);
      convergence = abs(f_x1);
      iters++;
    }
    if (iters == MAX_ITER) {
      throw nr_except;
    }
    result = x2;
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    result.resize(n, nan(""));
  }
  return result;
}

std::vector<double> constrained_solver_secant(
    std::vector<double> x0, double TOL,
    std::function<double(std::vector<double>)> objective,
    std::function<std::vector<double>(std::vector<double>)> constraints,
    const int MAX_ITER) {

  // allocate space for result and fill with zeros
  std::vector<double> result{};
  std::vector<double> x1 = x0;
  std::vector<double> x2;
  std::vector<double> lambda0{};
  std::vector<double> lambda1{};
  std::vector<double> lambda2;
  size_t n = x0.size();
  size_t c = constraints(x0).size();
  result.resize(n);
  for (int i = 0; i < n; i++) {
    x1[i] += TOL;
  }
  x2.resize(n);
  lambda0.resize(c, 1.0);
  lambda1.resize(c, 1.0 + TOL);
  lambda2.resize(c);
  // define the lagrangian of the system
  std::function<double(std::vector<double>, std::vector<double>)> L =
      [=](std::vector<double> x, std::vector<double> lambda) -> double {
    std::vector<double> g = constraints(x);
    return objective(x) -
           std::inner_product(g.begin(), g.end(), lambda.begin(), 0.0);
  };
  // iterate until convergence
  try {
    double convergence = 2 * TOL;
    int iters = 0;
    double L_x0 = L(x0, lambda0);
    double L_x1 = L(x1, lambda1);
    while (iters < MAX_ITER && TOL < convergence) {
      double step = L_x1 / (L_x1 - L_x0);
      for (int i = 0; i < n; i++) {
        x2[i] = x1[i] - step * (x1[i] - x0[i]);
      }
      for (int i = 0; i < c; i++) {
        lambda2[i] = lambda1[i] - step * (lambda1[i] - lambda0[i]);
      }
      x0 = x1;
      x1 = x2;
      lambda0 = lambda1;
      lambda1 = lambda2;
      L_x0 = L_x1;
      L_x1 = L(x1, lambda1);
      convergence = abs(L_x1);
      iters++;
    }
    if (iters == MAX_ITER) {
      throw nr_except;
    }
    result = x2;
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    result.resize(n, nan(""));
  }
  return result;
}

double RK4(std::function<double(double, double)> F, double yi, double xi,
           double dx) {
  double k1 = F(xi, yi);
  double k2 = F(xi + dx / 2.0, yi + k1 * dx / 2.0);
  double k3 = F(xi + dx / 2.0, yi + k2 * dx / 2.0);
  double k4 = F(xi + dx, yi + k3 * dx);
  return yi + (k1 + 2.0 * k2 + 2.0 * k3 + k4) * dx / 6.0;
}

double integrate(std::function<double(std::function<double(double, double)>,
                                      double, double, double)>
                     integration_method,
                 std::function<double(double, double)> F, double y_init,
                 double x_lower, double x_upper, double dx) {
  double y = y_init;
  double x = x_lower;
  while (x + dx <= x_upper) {
    y = integration_method(F, y, x, dx);
    x += dx;
  }
  if (x < x_upper) {
    dx = x_upper - x;
    y = integration_method(F, y, x, dx);
    x += dx;
  }
  return y;
}

double interp_1D(std::vector<std::pair<double, double>> &func, double x) {
  std::pair<double, double> p1;
  std::pair<double, double> p2;
  int last = func.size() - 1;
  assert(last >= 1);
  if (x <= func[0].first) {
    p1 = func[0];
    p2 = func[1];
  } else if (x >= func[last].first) {
    p1 = func[last - 1];
    p2 = func[last];
  } else {
    for (int i = 0; i <= last; i++) {
      if (func[i].first >= x) {
        p1 = func[i - 1];
        p2 = func[i];
        break;
      }
    }
  }
  double m = (p2.second - p1.second) / (p2.first - p1.first);
  double b = p1.second - m * p1.first;
  return m * x + b;
}

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

} // namespace hazen