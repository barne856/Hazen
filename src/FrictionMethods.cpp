#include "FrictionMethods.hpp"
#include "HydraulicUtil.hpp"
#include <assert.h>

namespace hazen {

ManningsFriction::ManningsFriction(double n) : n(n) { assert(n >= 0.0); }
double ManningsFriction::friction_slope(HydraulicShape *shape, double Q,
                                        double depth) {
  // Manning's Equation US units
  // Q = (k/n)*A*Rh^(2/3)*S^(1/2) --> S = (Q*n/(k*A*Rh^(2/3)))^2
  if (depth <= 0.0) {
    // return 0, positive, or negative infinity depending on flow direction
    if (Q > 0.0) {
      return std::numeric_limits<double>::infinity();
    }
    if (Q < 0.0) {
      return -std::numeric_limits<double>::infinity();
    }
    return 0.0;
  }
  if (depth == std::numeric_limits<double>::infinity()) {
    return 0.0;
  }

  const double k = 1.49;
  double A = shape->flow_area(depth);
  double Rh = shape->hydraulic_radius(depth);
  return Q * abs(Q) * pow((n / (k * A * pow(Rh, 2.0 / 3.0))), 2.0);
}

HazenFriction::HazenFriction(double C) : C(C) { assert(C >= 0.0); }
double HazenFriction::friction_slope(HydraulicShape *shape, double Q,
                                     double depth) {
  // Hazen-Williams Equation US units
  // Q = k*A*C*R^(0.63)*S^(0.54) --> S = (Q/(k*A*C*R^(0.63)))^(1/0.54)
  if (depth <= 0.0) {
    // return 0, positive, or negative infinity depending on flow direction
    if (Q > 0.0) {
      return std::numeric_limits<double>::infinity();
    }
    if (Q < 0.0) {
      return -std::numeric_limits<double>::infinity();
    }
    return 0.0;
  }
  if (depth == std::numeric_limits<double>::infinity()) {
    return 0.0;
  }

  double k = 1.318;
  double A = shape->flow_area(depth);
  double Rh = shape->hydraulic_radius(depth);
  return Q * abs(Q) * pow((1.0 / (k * A * C * pow(Rh, 0.63))), 1.0 / 0.54);
}

DarcyFriction::DarcyFriction(double eps, double rho, double mu)
    : eps(eps), rho(rho), mu(mu) {
  assert(eps >= 0.0);
  assert(rho >= 0.0);
  assert(mu >= 0.0);
}
double DarcyFriction::friction_slope(HydraulicShape *shape, double Q,
                                     double depth) {
  // Darcy-Weisbach Equation
  // S = f*V^2/(2g*Dh)
  if (depth <= 0.0) {
    // return 0, positive, or negative infinity depending on flow direction
    if (Q > 0.0) {
      return std::numeric_limits<double>::infinity();
    }
    if (Q < 0.0) {
      return -std::numeric_limits<double>::infinity();
    }
    return 0.0;
  }
  if (depth == std::numeric_limits<double>::infinity()) {
    return 0.0;
  }

  double A = shape->flow_area(depth);
  double Dh = shape->hydraulic_diameter(depth);
  double V = Q / A;
  double Re = reynolds(rho, mu, V, Dh);
  double f = shape->is_free_surface(depth)
                 ? darcy_friction_factor_free_surface(Re, Dh, eps)
                 : darcy_friction_factor_pressure_driven(Re, Dh, eps);
  return f * V * abs(V) / (2 * g * Dh);
}

} // namespace hazen