#include "Hazen/Hydraulics.hpp"
#include <assert.h>

namespace hazen {

ManningsFriction::ManningsFriction(Dimensionless n) : n(n.val) {
  assert(n.val >= 0.0);
}
Angle ManningsFriction::friction_slope(HydraulicShape *shape, Flow Q,
                                       Length depth) {
  // Manning's Equation US units
  // Q = (k/n)*A*Rh^(2/3)*S^(1/2) --> S = (Q*n/(k*A*Rh^(2/3)))^2
  const auto k = Scalar<Unit<1, 3, -1, 1, 0, 1, 0, 1>>(1.0);
  Area A = shape->flow_area(depth);
  Length Rh = shape->hydraulic_radius(depth);
  return Angle::Slope((Q * abs(Q) * abs2(n / (k * A * power<2, 3>(Rh)))).val,
                      1.0);
}

HazenFriction::HazenFriction(Dimensionless C) : C(C.val) {
  assert(C.val >= 0.0);
}
Angle HazenFriction::friction_slope(HydraulicShape *shape, Flow Q,
                                    Length depth) {
  // Hazen-Williams Equation US units
  // Q = k*A*C*R^(0.63)*S^(0.54) --> S = (Q/(k*A*C*R^(0.63)))^(1/0.54)
  const auto k = Scalar<Unit<37, 100, -1, 1, 0, 1, 0, 1>>(0.849);
  Area A = shape->flow_area(depth);
  Length Rh = shape->hydraulic_radius(depth);
  return Angle::Slope(
      (power<50, 27>((abs(Q) / (k * A * C * power<63, 100>(Rh))))).val, 1.0);
}

DarcyFriction::DarcyFriction(Length eps, Density rho, Dynamic_Viscosity mu)
    : eps(eps.val), rho(rho.val), mu(mu.val) {
  assert(eps.val >= 0.0);
  assert(rho.val >= 0.0);
  assert(mu.val >= 0.0);
}
Angle DarcyFriction::friction_slope(HydraulicShape *shape, Flow Q,
                                    Length depth) {
  // Darcy-Weisbach Equation
  // S = f*V^2/(2g*Dh)
  Area A = shape->flow_area(depth);
  Length Dh = shape->hydraulic_diameter(depth);
  Velocity V = Q / A;
  Dimensionless Re = reynolds(rho, mu, V, Dh);
  Dimensionless f = shape->is_free_surface(depth)
                        ? darcy_friction_factor_free_surface(Re, Dh, eps)
                        : darcy_friction_factor_pressure_driven(Re, Dh, eps);
  return Angle::Slope((f * V * abs(V) / (Dimensionless(2.0) * g * Dh)).val,
                      1.0);
}

} // namespace hazen