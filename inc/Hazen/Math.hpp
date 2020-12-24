#ifndef HAZEN_MATH
#define HAZEN_MATH

#include "Hazen/Core.hpp"

#include <functional>

namespace hazen {

/**
 * @brief Compute the horizontal length between two points.
 *
 * @param down downstream point
 * @param up upstream point
 * @return double The horizontal length between the two points.
 */
template <typename T>
T horizontal_length(const Vec<T> &down, const Vec<T> &up) {
  assert(down.size() == 3 && up.size() == 3);
  Vec<T> l(3);
  l.elems(0) = down(0).val - up(0).val;
  l.elems(1) = down(1).val - up(1).val;
  l.elems(2) = 0.0;
  return norm(l);
}

/**
 * @brief Compute the slope between and upstream and downstream point.
 * @details The slope is positive when the dowsntream point is lower than the
 * upstream point.
 *
 * @param down The downstream point
 * @param up The upstream point
 * @return double The slope from downstream to upstream.
 */
template <typename T> Angle slope(const Vec<T> &down, const Vec<T> &up) {
  assert(down.size() == 3 && up.size() == 3);
  T rise = up(2) - down(2);
  T run = horizontal_length(down, up);
  return Angle::Slope(rise.val, run.val);
}

// sum of squares for use in standard algorithms
template <typename T> struct sum_square {
  product_type<T, T> operator()(const product_type<T, T> &Left,
                                const T &Right) const {
    return (Left + Right * Right);
  }
};

// Root Finding Methods --------------------------------------------------------
/**
 * @brief Find the goal of an objective function using the Secant root
 * finding method.
 *
 * @param goal The goal of the objective function.
 * @param x0 The first initial starting point for the objective function. Should
 * ideally be chosen to lie close to the goal.
 * @param x1 The second initial starting point for the objective function.
 * Should ideally be chosen to lie close to the goal.
 * @param TOL The tolerance of the convergence for the solution.
 * @param objective The objective function.
 * @param MAX_ITER The maximum iterations to perform before the function quits
 * trying to refine the estimate.
 * @return double The input to the objective function that results in an
 * approximate value of goal being returned. NaN is returned if the goal cannot
 * be found or if the maximum iterations is exceeded.
 */
template <typename T, typename U>
T find_goal_secant(U goal, T x0, T x1, double TOL,
                   std::function<U(T)> objective, const int MAX_ITER = 1000) {
  T result{};

  double convergence = 2 * TOL;
  int iters = 0;
  T x2{};
  U f_x0 = objective(x0) - goal;
  U f_x1 = objective(x1) - goal;
  while (iters < MAX_ITER && TOL < convergence) {
    x2 = x1 - f_x1 * (x1 - x0) / (f_x1 - f_x0);
    x0 = x1;
    x1 = x2;
    f_x0 = f_x1;
    f_x1 = objective(x1) - goal;
    convergence = std::abs(f_x1.val);
    iters++;
  }
  result = x2;
  if (iters == MAX_ITER) {
    result = T(nan(""));
  }
  return result;
}

/**
 * @brief Find the goal of an objective function using the Bisection root
 * finding method.
 *
 * @param goal The goal of the objective function.
 * @param a The first initial starting point for the objective function. Should
 * be chosen with "b" to bracket to root. objective(a)*objective(b) must be
 * less than 0.
 * @param b The second initial starting point for the objective function.
 * Should be chosen with "a" to bracket to root. objective(a)*objective(b) must
 * be less than 0.
 * @param TOL The tolerance of the convergence for the solution.
 * @param objective The objective function.
 * @param MAX_ITER The maximum iterations to perform before the function quits
 * trying to refine the estimate.
 * @return double The input to the objective function that results in an
 * approximate value of goal being returned. NaN is returned if the goal cannot
 * be found or if the maximum iterations is exceeded.
 */
template <typename T, typename U>
T find_goal_bisection(U goal, T a, T b, double TOL,
                      std::function<U(T)> objective,
                      const int MAX_ITER = 1000) {
  T result{};
  int i = 0;
  T c{};
  U fc{};
  while (i < MAX_ITER) {
    c = (a + b) / 2.0;
    fc = objective(c) - goal;
    if (fc.val == 0.0 || (b - c).val / 2.0 < TOL) {
      // solution found
      break;
    }
    (fc * (objective(a) - goal)).val > 0.0 ? a = c : b = c;
    i++;
  }
  result = c;
  if (i == MAX_ITER) {
    result = T(nan(""));
  }
  return result;
}

/**
 * @brief Find the goal of an objective function using the Inverse Quadratic
 * root finding method.
 *
 * @param goal The goal of the objective function.
 * @param x0 The first initial starting point for the objective function. Should
 * ideally be chosen to lie close to the goal.
 * @param x1 The second initial starting point for the objective function.
 * Should ideally be chosen to lie close to the goal.
 * @param x2 The third initial starting point for the objective function.
 * Should ideally be chosen to lie close to the goal.
 * @param TOL The tolerance of the convergence for the solution.
 * @param objective The objective function.
 * @param MAX_ITER The maximum iterations to perform before the function quits
 * trying to refine the estimate.
 * @return double The input to the objective function that results in an
 * approximate value of goal being returned. NaN is returned if the goal cannot
 * be found or if the maximum iterations is exceeded.
 */
template <typename T, typename U>
Scalar<T> find_goal_inverse_quadratic(
    Scalar<U> goal, Scalar<T> x0, Scalar<T> x1, Scalar<T> x2, double TOL,
    std::function<Scalar<U>(Scalar<T>)> objective, const int MAX_ITER = 1000) {
  Scalar<T> result{};
  double convergence = 2 * TOL;
  int iters = 0;
  Scalar<T> x3{};
  Scalar<U> f_x0 = objective(x0) - goal;
  Scalar<U> f_x1 = objective(x1) - goal;
  Scalar<U> f_x2 = objective(x2) - goal;
  while (iters < MAX_ITER && TOL < convergence) {
    Dimensionless a = f_x1 * f_x2 / ((f_x0 - f_x1) * (f_x0 - f_x2));
    Dimensionless b = f_x0 * f_x2 / ((f_x1 - f_x0) * (f_x1 - f_x2));
    Dimensionless c = f_x0 * f_x1 / ((f_x2 - f_x0) * (f_x2 - f_x1));
    x3 = a * x0 + b * x1 + c * x2;
    x0 = x1;
    x1 = x2;
    x2 = x3;
    f_x0 = f_x1;
    f_x1 = f_x2;
    f_x2 = objective(x2) - goal;
    convergence = abs(f_x2.val);
    iters++;
  }
  result = x3;
  if (iters == MAX_ITER) {
    result = Scalar<T>(nan(""));
  }
  return result;
}
// Linear Algebra --------------------------------------------------------------
template <typename T, typename U>
Vec<quotient_type<U, T>> gradient(Vec<T> x0, std::function<U(Vec<T>)> objective,
                                  T dx) {
  Vec<quotient_type<U, T>> result(x0.size());
  U f0 = objective(x0);
  for (int i = 0; i < x0.size(); i++) {
    Vec<T> x1 = x0;
    x1(i) = x1(i) + dx;
    U f1 = objective(x1);
    result(i) = (f1 - f0) / dx;
  }
  return result;
}

template <typename T, typename U>
Mat<quotient_type<U, T>>
Jacobian(Vec<T> x0, std::function<Vec<U>(Vec<T>)> objective, T dx) {
  Vec<U> F0 = objective(x0);
  int n = x0.size();
  int m = F0.size();
  Mat<quotient_type<U, T>> JF(m, n);
  for (int j = 0; j < n; j++) {
    Vec<T> x_new = x0;
    x_new(j) = x_new(j) + dx;
    Vec<U> Fj = objective(x_new);
    for (int i = 0; i < m; i++) {
      JF(i, j) = (Fj(i) - F0(i)) / dx;
    }
  }
  return JF;
}

template <typename T, typename U>
Vec<T>
solve_least_squares_linear_system(Vec<T> x0, double TOL,
                                  std::function<Vec<U>(Vec<T>)> objective) {
  // compute objective at x0
  Vec<U> F_x0 = objective(x0);
  // compute Jacobian at x0
  auto JF = Jacobian(x0, objective, T(TOL));
  // compute delta x using the moore-penrose pseudoinverse
  Vec<T> delta_x = pinv(JF) * F_x0;
  return x0 - delta_x;
}

template <typename T, typename U>
Vec<T> find_root_generalized_newton(Vec<T> x0,
                                    std::function<Vec<U>(Vec<T>)> objective,
                                    double TOL, const int MAX_ITER = 1000) {
  // iterate until convergence
  double convergence = 2 * TOL;
  int iters = 0;
  Vec<U> F_x0 = objective(x0);
  while (iters < MAX_ITER && TOL < convergence) {
    auto JF = Jacobian(x0, objective, T(TOL));
    x0 -= solve(JF, F_x0);
    F_x0 = objective(x0);
    convergence = std::accumulate(F_x0.begin(), F_x0.end(),
                                  product_type<U, U>(0.0), sum_square<U>())
                      .val;
    iters++;
  }
  if (iters == MAX_ITER) {
    x0 = T(nan(""));
  }
  return x0;
}

template <typename T, typename U, typename V>
Vec<T>
solve_constrained_problem_lagrange(Vec<T> x0,
                                   std::function<U(Vec<T>)> objective,
                                   std::function<Vec<V>(Vec<T>)> constraints,
                                   double TOL, const int MAX_ITER = 1000) {
  std::function<U(Vec<T>)> lagrangian = [=](Vec<T> s) {
    auto &x = s(slice(0, x0.size()));
    auto &lambda = s(slice(x0.size()));
    U f = objective(x);
    Vec<V> g = constraints(x);
    U result = f - dot_product(lambda, g(slice(0))) *
                       quotient_type<U, product_type<T, V>>(1.0);
    return result;
  };
  std::function<Vec<quotient_type<U, T>>(Vec<T>)> grad_lagrangian =
      [=](Vec<T> s) { return gradient(s, lagrangian, T(TOL)); };
  Vec<T> c = constraints(x0) * quotient_type<T, V>(1.0);
  Vec<T> s(x0.size() + c.size());
  s(slice(0, x0.size())) = x0;
  s(slice(x0.size())) = c;
  auto solution =
      find_root_generalized_newton(s, grad_lagrangian, TOL, MAX_ITER);
  return solution(slice(0, x0.size()));
}

// Numerical Integration -------------------------------------------------------

namespace RKF45_CONSTANTS {
const Dimensionless a[6] = {
    Dimensionless(0.0),        Dimensionless(1.0 / 5.0),
    Dimensionless(3.0 / 10.0), Dimensionless(3.0 / 5.0),
    Dimensionless(1.0),        Dimensionless(7.0 / 8.0)};
const Dimensionless b[6][5] = {
    {Dimensionless(0.0), Dimensionless(0.0), Dimensionless(0.0),
     Dimensionless(0.0), Dimensionless(0.0)},
    {Dimensionless(1.0 / 5.0), Dimensionless(0.0), Dimensionless(0.0),
     Dimensionless(0.0), Dimensionless(0.0)},
    {Dimensionless(3.0 / 40.0), Dimensionless(9.0 / 40.0), Dimensionless(0.0),
     Dimensionless(0.0), Dimensionless(0.0)},
    {Dimensionless(3.0 / 10.0), Dimensionless(-9.0 / 10.0),
     Dimensionless(6.0 / 5.0), Dimensionless(0.0), Dimensionless(0.0)},
    {Dimensionless(-11.0 / 54.0), Dimensionless(5.0 / 2.0),
     Dimensionless(-70.0 / 27.0), Dimensionless(35.0 / 27.0),
     Dimensionless(0.0)},
    {Dimensionless(1631.0 / 55296.0), Dimensionless(175.0 / 512.0),
     Dimensionless(575.0 / 13824.0), Dimensionless(44275.0 / 110592.0),
     Dimensionless(253.0 / 4096.0)}};
const Dimensionless c[6] = {
    Dimensionless(37.0 / 378.0),  Dimensionless(0.0),
    Dimensionless(250.0 / 621.0), Dimensionless(125.0 / 594.0),
    Dimensionless(0.0),           Dimensionless(512.0 / 1771.0)};
const Dimensionless d[6] = {
    Dimensionless(c[0].val) - Dimensionless(2825.0 / 27648.0),
    Dimensionless(c[1].val) - Dimensionless(0.0),
    Dimensionless(c[2].val) - Dimensionless(18575.0 / 48384.0),
    Dimensionless(c[3].val) - Dimensionless(13525.0 / 55296.0),
    Dimensionless(c[4].val) - Dimensionless(277.0 / 14336.0),
    Dimensionless(c[5].val) - Dimensionless(1.0 / 4.0)};
} // namespace RKF45_CONSTANTS

/**
 * @brief Compute the solution of an ODE using a Runge-Kutta (4,5) embedded
 * algorithm with an adaptive step size.
 *
 * @param ode The ODE to solve dy/dx(x,y).
 * @param init The initial value y0.
 * @param span The span [x0, xf] to solve for.
 * @param rel_tol The relative error tolerance criterion used to compute the
 * adaptive step size.
 * @return std::vector<double> The result y(x) on the span [x0, xf].
 */
template <typename X, typename Y>
std::vector<std::pair<X, Y>>
RKF45(std::function<quotient_type<Y, X>(X, Y)> ode, Y init,
      std::pair<X, X> span, std::function<bool(X, Y)> condition = nullptr,
      Y rel_tol = Y(1e-5), int MAX_STEP = 1000, X h_min = X(1e-5)) {
  using namespace RKF45_CONSTANTS;
  // relative tolerance must by strictly positive real number.
  assert(rel_tol > Y(0.0) && !isinf(rel_tol.val) && !isnan(rel_tol.val));
  // initialize variables
  std::vector<std::pair<X, Y>> result{};
  X h = span.second - span.first;
  X x = span.first;
  Y y = init;
  result.push_back({x, y});
  double sign = h.val > 0.0 ? 1.0 : -1.0;
  int step = 0;
  // loop until the end of the span.
  while (sign * x.val < sign * span.second.val) {
    step += 1;
    if (step > MAX_STEP) {
      h = X(nan(""));
    }
    if (abs(h) < h_min) {
      h = X(nan(""));
    }
    if (isnan(h.val)) {
      break;
    }
    if (sign * (x.val + h.val) > sign * span.second.val + rel_tol.val) {
      // result will be beyond the end of the span, set to the end.
      h = span.second - x;
      continue;
    }
    // RK fifth order k values
    Y k[6] = {h * ode(x, y),
              h * ode(x + a[1] * h, y + b[1][0] * k[0]),
              h * ode(x + a[2] * h, y + b[2][0] * k[0] + b[2][1] * k[1]),
              h * ode(x + a[3] * h,
                      y + b[3][0] * k[0] + b[3][1] * k[1] + b[3][2] * k[2]),
              h * ode(x + a[4] * h, y + b[4][0] * k[0] + b[4][1] * k[1] +
                                        b[4][2] * k[2] + b[4][3] * k[3]),
              h * ode(x + a[5] * h, y + b[5][0] * k[0] + b[5][1] * k[1] +
                                        b[5][2] * k[2] + b[5][3] * k[3] +
                                        b[5][4] * k[4])};
    // relative error between the 4th and 5th order methods.
    Y delta =
        d[0] * k[0] + d[2] * k[2] + d[3] * k[3] + d[4] * k[4] + d[5] * k[5];
    // adjust step size and record values if they are valid.
    if (!isnan(delta.val) && abs(delta) <= abs(rel_tol)) {
      x += h;
      y += c[0] * k[0] + c[2] * k[2] + c[3] * k[3] + c[5] * k[5];
      result.push_back({x, y});
      if (condition && condition(x, y) == false) {
        result.pop_back();
        h = X(nan(""));
      }
      // if delta is zero, the solution is exact and we can skip to the end.
      if (delta.val == 0.0) {
        h = span.second - x;
      } else {
        h = X(0.95 * h.val * std::pow(std::abs(rel_tol.val / delta.val), 0.2));
      }
    } else {
      h = X(0.95 * h.val * std::pow(std::abs(rel_tol.val / delta.val), 0.25));
    }
  }
  if (isnan(h.val)) {
    result.push_back({x, Y(nan(""))});
  }
  return result;
}

/**
 * @brief Linearly interpolate a function of a single variable given as a vector
 * of x,y pairs.
 *
 * @param data The data to interpolate, vector of (x,y) pairs.
 * @param x The x value to interpolate
 * @return double The interpolated y value at x
 */
template <typename U1, typename U2>
U2 interp_1D(std::vector<std::pair<U1, U2>> &data, U1 x) {
  std::pair<U1, U2> p1;
  std::pair<U1, U2> p2;
  int last = data.size() - 1;
  assert(last >= 1);
  if (x <= data[0].first) {
    p1 = data[0];
    p2 = data[1];
  } else if (x >= data[last].first) {
    p1 = data[last - 1];
    p2 = data[last];
  } else {
    for (int i = 0; i <= last; i++) {
      if (data[i].first >= x) {
        p1 = data[i - 1];
        p2 = data[i];
        break;
      }
    }
  }
  double m = (p2.second - p1.second) / (p2.first - p1.first);
  double b = p1.second - m * p1.first;
  return m * x + b;
}

} // namespace hazen

#endif