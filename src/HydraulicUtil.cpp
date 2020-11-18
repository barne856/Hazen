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

template <typename T> struct square {
  T operator()(const T &Left, const T &Right) const {
    return (Left + Right * Right);
  }
};

std::vector<std::vector<double>>
Jacobian(std::vector<double> x0, std::vector<double> F0,
         std::function<std::vector<double>(std::vector<double>)> objective,
         double TOL) {
  int n = x0.size();
  int m = F0.size();
  std::vector<std::vector<double>> JF{};
  JF.resize(m);
  for (int i = 0; i < m; i++) {
    JF[i].resize(n);
  }
  for (int j = 0; j < n; j++) {
    std::vector<double> x_new = x0;
    x_new[j] += TOL;
    std::vector<double> Fj = objective(x_new);
    for (int i = 0; i < m; i++) {
      JF[i][j] = (Fj[i] - F0[i]) / TOL;
    }
  }
  return JF;
}

std::vector<double>
gradient(std::vector<double> x0,
         std::function<double(std::vector<double>)> objective, double TOL) {
  std::vector<double> result{};
  result.reserve(x0.size());
  double f0 = objective(x0);
  for (int i = 0; i < x0.size(); i++) {
    auto x1 = x0;
    x1[i] += TOL;
    double f1 = objective(x1);
    result.push_back((f1 - f0) / TOL);
  }
  return result;
}

typedef std::vector<std::vector<double>> mat;
typedef std::vector<double> vec;

void print_vec(vec vector) {
  for (auto &x : vector) {
    std::cout << x << ", ";
  }
  std::cout << std::endl;
}

void print_mat(mat matrix) {
  for (auto &vector : matrix) {
    print_vec(vector);
  }
}

void pivot(mat &Ab) {
  int m = Ab.size();
  int n = Ab[0].size() - 1;
  for (int i = m - 1; i > 0; i--) // partial pivoting
  {
    if (Ab[i - 1][0] < Ab[i][0])
      for (int j = 0; j <= n; j++) {
        double c = Ab[i][j];
        Ab[i][j] = Ab[i - 1][j];
        Ab[i - 1][j] = c;
      }
  }
}

void forward_elim(mat &Ab) {
  int m = Ab.size();
  int n = Ab[0].size() - 1;
  for (int k = 0; k < m - 1; k++)
    for (int i = k; i < m - 1; i++) {
      double c = (Ab[i + 1][k] / Ab[k][k]);

      for (int j = 0; j <= n; j++)
        Ab[i + 1][j] -= c * Ab[k][j];
    }
}

vec back_sub(const mat &Ab) {
  vec x{}; // A vec to store solution
  int m = Ab.size();
  int n = Ab[0].size() - 1;
  x.resize(m);
  for (int i = m - 1; i >= 0; i--) {
    double c = 0;
    for (int j = i; j <= n - 1; j++)
      c = c + Ab[i][j] * x[j];

    x[i] = (Ab[i][n] - c) / Ab[i][i];
  }
  return x;
}

vec solve_linear_system(const mat &A, const vec &b) {
  // create augmented matrix Ab
  mat Ab{};
  for (int i = 0; i < A.size(); i++) {
    vec row = A[i];
    row.push_back(b[i]);
    Ab.push_back(row);
  }
  // partial pivoting
  pivot(Ab);
  // reduce to r.e.f.
  forward_elim(Ab);
  // solve for x
  return back_sub(Ab);
}

vec gauss(const mat &A, const vec &b) {
  mat Ab{};
  for (int i = 0; i < A.size(); i++) {
    vec row = A[i];
    row.push_back(b[i]);
    Ab.push_back(row);
  }
  int m = Ab.size();
  int n = Ab[0].size() - 1;

  for (int i = 0; i < m; i++) {
    // Search for maximum in this column
    double maxEl = abs(Ab[i][i]);
    int maxRow = i;
    for (int k = i + 1; k < m; k++) {
      if (abs(Ab[k][i]) > maxEl) {
        maxEl = abs(Ab[k][i]);
        maxRow = k;
      }
    }

    // Swap maximum row with current row (column by column)
    for (int k = i; k < n + 1; k++) {
      double tmp = Ab[maxRow][k];
      Ab[maxRow][k] = Ab[i][k];
      Ab[i][k] = tmp;
    }

    // Make all rows below this one 0 in current column
    for (int k = i + 1; k < m; k++) {
      double c = -Ab[k][i] / Ab[i][i];
      for (int j = i; j < n + 1; j++) {
        if (i == j) {
          Ab[k][j] = 0;
        } else {
          Ab[k][j] += c * Ab[i][j];
        }
      }
    }
  }

  // Solve equation Ax=b for an upper triangular matrix A
  vec x(n);
  for (int i = m - 1; i >= 0; i--) {
    x[i] = Ab[i][n] / Ab[i][i];
    for (int k = i - 1; k >= 0; k--) {
      Ab[k][n] -= Ab[k][i] * x[i];
    }
  }
  return x;
}

mat transpose(const mat &A) {
  int m = A.size();
  int n = A[0].size();
  // allocate space for B = A^T
  mat B(n);
  for (int i = 0; i < n; i++) {
    B[i].resize(m);
  }
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < m; j++) {
      B[i][j] = A[j][i];
    }
  }
  return B;
}

mat multiply(const mat &A, const mat &B) {
  int m = A.size();
  int n = A[0].size();
  int o = B[0].size();
  // allocate space for C = A x B
  mat C(m);
  for (int i = 0; i < m; i++) {
    C[i].resize(o);
  }
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < o; j++) {
      double d = 0.0;
      for (int k = 0; k < n; k++) {
        d += A[i][k] * B[k][j];
      }
      C[i][j] = d;
    }
  }
  return C;
}

mat inverse(const mat &A) {
  // create augmented matrix AB
  mat AB{};
  int n = A.size();
  for (int i = 0; i < n; i++) {
    vec row = A[i];
    for (int j = 0; j < n; j++) {
      if (i == j) {
        row.push_back(1.0);
      } else {
        row.push_back(0.0);
      }
    }
    AB.push_back(row);
  }

  for (int i = 0; i < n; i++) {
    // Search for maximum in this column
    double maxEl = abs(AB[i][i]);
    int maxRow = i;
    for (int k = i + 1; k < n; k++) {
      if (abs(AB[k][i]) > maxEl) {
        maxEl = abs(AB[k][i]);
        maxRow = k;
      }
    }

    // Swap maximum row with current row (column by column)
    for (int k = i; k < 2 * n; k++) {
      double tmp = AB[maxRow][k];
      AB[maxRow][k] = AB[i][k];
      AB[i][k] = tmp;
    }

    // Make all rows below this one 0 in current column
    for (int k = i + 1; k < n; k++) {
      double c = -AB[k][i] / AB[i][i];
      for (int j = i; j < 2 * n; j++) {
        if (i == j) {
          AB[k][j] = 0;
        } else {
          AB[k][j] += c * AB[i][j];
        }
      }
    }
  }

  // Solve equation AB=I for an upper triangular matrix A
  mat B(n);
  for (int i = 0; i < n; i++) {
    B[i].resize(n);
  }
  for (int i = n - 1; i >= 0; i--) {
    for (int j = 0; j < n; j++) {
      B[i][j] = AB[i][n + j] / AB[i][i];
    }
    for (int k = i - 1; k >= 0; k--) {
      for (int j = 0; j < n; j++) {
        AB[k][n + j] -= AB[k][i] * B[i][j];
      }
    }
  }
  return B;
}

mat pseudoinverse(const mat &A) {
  mat A_T = transpose(A);
  return multiply(A_T, inverse(multiply(A, A_T)));
}

std::vector<double> solve_least_squares_linear_system(
    std::vector<double> x0, double TOL,
    std::function<std::vector<double>(std::vector<double>)> objective) {
  // allocate space for result and fill with zeros
  std::vector<double> result = x0;
  // compute objective at x0
  std::vector<double> F_x0 = objective(x0);
  // compute Jacobian at x0
  std::vector<std::vector<double>> JF = Jacobian(x0, F_x0, objective, TOL);
  // compute delta x using the moore-penrose pseudoinverse
  mat b{};
  b.push_back(F_x0);
  b = transpose(b);
  mat delta_x = multiply(pseudoinverse(JF), b);
  // adjust result
  for (int i = 0; i < x0.size(); i++) {
    result[i] -= delta_x[i][0];
  }
  return result;
}

std::vector<double> find_root_generalized_newton(
    std::vector<double> x0,
    std::function<std::vector<double>(std::vector<double>)> objective,
    double TOL, const int MAX_ITER) {
  // iterate until convergence
  try {
    double convergence = 2 * TOL;
    int iters = 0;
    std::vector<double> F_x0 = objective(x0);
    while (iters < MAX_ITER && TOL < convergence) {
      std::vector<std::vector<double>> JF = Jacobian(x0, F_x0, objective, TOL);
      vec delta_x = gauss(JF, F_x0);
      for (int i = 0; i < x0.size(); i++) {
        x0[i] -= delta_x[i];
      }
      F_x0 = objective(x0);
      convergence =
          std::accumulate(F_x0.begin(), F_x0.end(), 0.0, square<double>());
      iters++;
    }
    if (iters == MAX_ITER) {
      throw nr_except;
    }
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    x0.resize(x0.size(), nan(""));
  }
  return x0;
}

std::vector<double> solve_constrained_problem_lagrange(
    std::vector<double> x0,
    std::function<double(std::vector<double>)> objective,
    std::function<std::vector<double>(std::vector<double>)> constraints,
    double TOL, const int MAX_ITER) {
  std::function<double(std::vector<double>)> lagrangian =
      [=](std::vector<double> s) {
        auto &x = std::vector<double>(s.begin(), s.begin() + x0.size());
        auto &lambda = std::vector<double>(s.begin() + x0.size(), s.end());
        double f = objective(x);
        std::vector<double> g = constraints(x);
        double result = f;
        for (int i = 0; i < lambda.size(); i++) {
          result -= lambda[i] * g[i];
        }
        return result;
      };
  std::function<std::vector<double>(std::vector<double>)> grad_lagrangian =
      [=](std::vector<double> s) { return gradient(s, lagrangian, TOL); };
  auto c = constraints(x0);
  c.resize(c.size(), 1.0);
  auto s = x0;
  s.insert(s.end(), c.begin(), c.end());
  auto solution =
      find_root_generalized_newton(s, grad_lagrangian, TOL, MAX_ITER);
  return std::vector<double>(solution.begin(), solution.begin() + x0.size());
}

std::vector<double> secant_root_multivariate(
    std::vector<double> x0, double TOL,
    std::function<std::vector<double>(std::vector<double>)> objective,
    const int MAX_ITER) {
  // allocate space for result and fill with zeros
  size_t n = x0.size();
  std::vector<double> result{};
  result.resize(n);
  // allocate space for x1 and initialize to zero
  std::vector<double> x1 = x0;

  // iterate until convergence
  try {
    double convergence = 2 * TOL;
    int iters = 0;
    std::vector<double> F_x0 = objective(x0);
    while (iters < MAX_ITER && TOL < convergence) {
      std::vector<std::vector<double>> JF = Jacobian(x0, F_x0, objective, TOL);
      mat b{};
      b.push_back(F_x0);
      b = transpose(b);
      mat delta_x = multiply(pseudoinverse(JF), b);
      for (int i = 0; i < x1.size(); i++) {
        x1[i] -= delta_x[i][0];
      }
      x0 = x1;
      F_x0 = objective(x0);
      convergence =
          std::accumulate(F_x0.begin(), F_x0.end(), 0.0, square<double>());
      iters++;
    }
    if (iters == MAX_ITER) {
      throw nr_except;
    }
    result = x1;
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

namespace RKF45_CONSTANTS {
const double a[6] = {0.0, 1.0 / 5.0, 3.0 / 10.0, 3.0 / 5.0, 1.0, 7.0 / 8.0};
const double b[6][5] = {
    {0.0, 0.0, 0.0, 0.0, 0.0},
    {1.0 / 5.0, 0.0, 0.0, 0.0, 0.0},
    {3.0 / 40.0, 9.0 / 40.0, 0.0, 0.0, 0.0},
    {3.0 / 10.0, -9.0 / 10.0, 6.0 / 5.0, 0.0, 0.0},
    {-11.0 / 54.0, 5.0 / 2.0, -70.0 / 27.0, 35.0 / 27.0, 0.0},
    {1631.0 / 55296.0, 175.0 / 512.0, 575.0 / 13824.0, 44275.0 / 110592.0,
     253.0 / 4096.0}};
const double c[6] = {37.0 / 378.0,  0.0, 250.0 / 621.0,
                     125.0 / 594.0, 0.0, 512.0 / 1771.0};
const double d[6] = {c[0] - 2825.0 / 27648.0,  c[1] - 0.0,
                     c[2] - 18575.0 / 48384.0, c[3] - 13525.0 / 55296.0,
                     c[4] - 277.0 / 14336.0,   c[5] - 1.0 / 4.0};
} // namespace RKF45_CONSTANTS

std::vector<std::pair<double, double>>
RKF45(std::function<double(double, double)> ode, double init,
      std::pair<double, double> span,
      std::function<bool(double, double)> condition, double rel_tol,
      int MAX_STEP, double h_min) {
  using namespace RKF45_CONSTANTS;
  // relative tolerance must by strictly positive real number.
  assert(rel_tol > 0.0 && !isinf(rel_tol) && !isnan(rel_tol));
  // initialize variables
  std::vector<std::pair<double, double>> result{};
  double h = span.second - span.first;
  double x = span.first;
  double y = init;
  result.push_back({x, y});
  double sign = h > 0 ? 1.0 : -1.0;
  int step = 0;
  // loop until the end of the span.
  while (sign * x < sign * span.second) {
    step += 1;
    if (step > MAX_STEP) {
      h = nan("");
    }
    if (abs(h) < h_min) {
      h = nan("");
    }
    if (isnan(h)) {
      break;
    }
    if (sign * (x + h) > sign * span.second + rel_tol) {
      // result will be beyond the end of the span, set to the end.
      h = span.second - x;
      continue;
    }
    // RK fifth order k values
    double k[6] = {h * ode(x, y),
                   h * ode(x + a[1] * h, y + b[1][0] * k[0]),
                   h * ode(x + a[2] * h, y + b[2][0] * k[0] + b[2][1] * k[1]),
                   h * ode(x + a[3] * h, y + b[3][0] * k[0] + b[3][1] * k[1] +
                                             b[3][2] * k[2]),
                   h * ode(x + a[4] * h, y + b[4][0] * k[0] + b[4][1] * k[1] +
                                             b[4][2] * k[2] + b[4][3] * k[3]),
                   h * ode(x + a[5] * h, y + b[5][0] * k[0] + b[5][1] * k[1] +
                                             b[5][2] * k[2] + b[5][3] * k[3] +
                                             b[5][4] * k[4])};
    // relative error between the 4th and 5th order methods.
    double delta =
        d[0] * k[0] + d[2] * k[2] + d[3] * k[3] + d[4] * k[4] + d[5] * k[5];
    // adjust step size and record values if they are valid.
    if (!isnan(delta) && abs(delta) <= abs(rel_tol)) {
      x += h;
      y += c[0] * k[0] + c[2] * k[2] + c[3] * k[3] + c[5] * k[5];
      result.push_back({x, y});
      if (condition && condition(x, y) == false) {
        result.pop_back();
        h = nan("");
      }
      // if delta is zero, the solution is exact and we can skip to the end.
      if (delta == 0.0) {
        h = span.second - x;
      } else {
        h = 0.95 * h * pow(abs(rel_tol / delta), 0.2);
      }
    } else {
      h = 0.95 * h * pow(abs(rel_tol / delta), 0.25);
    }
  }
  if (isnan(h)) {
    result.push_back({x, nan("")});
  }
  return result;
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

std::vector<std::pair<std::string, std::vector<double>>>
gen_table(std::vector<std::pair<double, double>> values, std::string label1,
          std::string label2) {
  std::vector<double> vec1{};
  std::vector<double> vec2{};
  for (auto [x, y] : values) {
    vec1.push_back(x);
    vec2.push_back(y);
  }
  std::vector<std::pair<std::string, std::vector<double>>> cols{};
  cols.push_back({label1, vec1});
  cols.push_back({label2, vec2});
  return cols;
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