#include "HydraulicUtil.hpp"
#include <algorithm>
#include <assert.h>
#include <exception>
#include <iostream>
#include <math.h>

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

double darcy_friction_factor(double Re, double Dh, double eps) {
  // Explicit approximation valid for all flow regimes from Bellos, Nalbantis,
  // and Tsakiris (2018).
  double a = 1.0 / (1.0 + pow(Re / 2712.0, 8.4));
  double b = 1.0 / (1.0 + pow(Re / (150.0 * Dh / eps), 1.8));
  double f1 = 64.0 / Re;
  double f2 = 0.75 * log(Re / 5.37);
  double f3 = 0.88 * log(6.82 * Dh / eps);
  return pow(f1, a) * pow(f2, 2.0 * (a - 1.0) * b) *
         pow(f3, 2.0 * (a - 1.0) * (1.0 - b));
}

double normal_depth(HydraulicShape *shape, FrictionMethod *friction, double S,
                    double Q) {
  if (Q == 0.0) {
    return std::numeric_limits<double>::infinity();
  }
  if (abs(friction->friction_slope(shape, Q, shape->get_max_depth())) >=
      abs(S)) {
    return std::numeric_limits<double>::infinity();
  }
  // define function to use with bisection method
  std::function<double(double)> objective = [shape, friction,
                                             Q](double depth) -> double {
    return friction->friction_slope(shape, Q, depth);
  };
  double TOL = 0.00000001 * shape->height;
  double a = 0.5 * shape->height;
  double b = a;
  while (objective(a) - S < 0.0) {
    a /= 2.0;
    if (a < TOL) {
      return 0.0;
    }
  }
  while (objective(b) - S > 0.0) {
    b *= 2;
    if (b > shape->get_max_depth()) {
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
  double TOL = 0.00000001 * shape->height;
  double fa = (objective(TOL) - 1.0);
  double fb = (objective(shape->get_max_depth()) - 1.0);
  if (fa * fb > 0.0) {
    return std::numeric_limits<double>::infinity();
  }
  double a = TOL;
  double b = shape->get_max_depth();
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

double clamp(double value, double min, double max) {
  if (value < min) {
    return min;
  }
  if (value > max) {
    return max;
  }
  return value;
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

double alignment_length(alignment &align) {
  double L = 0.0;
  for (int i = 0; i < align.size() - 1; i++) {
    L += length(align[i], align[i + 1]);
  }
  return L;
}

double alignment_horizontal_length(alignment &align) {
  double L = 0.0;
  for (int i = 0; i < align.size() - 1; i++) {
    L += horizontal_length(align[i], align[i + 1]);
  }
  return L;
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
    double convergence = std::numeric_limits<double>::infinity();
    int i = 0;
    double x2;
    double fx0 = objective(x0) - goal;
    double fx1 = objective(x1) - goal;
    while (i < MAX_ITER && TOL < convergence) {
      x2 = x1 - fx1 * (x1 - x0) / (fx1 - fx0);
      x0 = x1;
      x1 = x2;
      fx0 = fx1;
      fx1 = objective(x1) - goal;
      convergence = abs(fx1);
      i++;
    }
    if (i == MAX_ITER) {
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
    double convergence = std::numeric_limits<double>::infinity();
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

double RK4(std::function<double(double, double)> F, double yi, double xi,
           double dx) {
  double k1 = F(xi, yi);
  double k2 = F(xi + dx / 2.0, yi + k1 * dx / 2.0);
  double k3 = F(xi + dx / 2.0, yi + k2 * dx / 2.0);
  double k4 = F(xi + dx, yi + k3 * dx);
  return yi + (k1 + 2.0 * k2 + 2.0 * k3 + k4) * dx / 6.0;
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

} // namespace hazen