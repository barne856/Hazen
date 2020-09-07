#ifndef FRICTIONMETHODS
#define FRICTIONMETHODS

// HAZEN
#include "HydraulicShapes.hpp"

namespace hazen {
/**
 * @brief An Abstract Friction Method used to calculate head losses in the
 * Conduit Hydraulic Link..
 */
class FrictionMethod {
public:
  /**
   * @brief Calculates the slope of the Energy Grade Line (EGL) at one point in
   * a Conduit.
   *
   * @param shape The cross-sectional shape of the Conduit.
   * @param Q The flow through the Conduit [UNITS = CFS].
   * @param depth The depth of flow in the Conduit at the calculation point
   * [UNITS = FT].
   * @return double The slope of the EGL at a point along the Conduit.
   */
  virtual double friction_slope(HydraulicShape *shape, double Q,
                                double depth) = 0;
};

/**
 * @brief An implementation of Manning's Equation to calculate the slope of the
 * Energy Grade Line (EGL) in a Conduit.
 *
 */
class ManningsFriction : public FrictionMethod {
public:
  /**
   * @brief Construct a new Mannings Friction object
   *
   * @param n Manning's n coefficient.
   */
  ManningsFriction(double n);
  /**
   * @brief An implementation of Manning's Equation to calculate the slope of
   * the Energy Grade Line (EGL) in a Conduit.
   * @param shape The cross-sectional shape of the Conduit.
   * @param Q The flow through the Conduit [UNITS = CFS].
   * @param depth The depth of flow in the Conduit at the calculation point
   * [UNITS = FT].
   * @return double The slope of the EGL at a point along the Conduit [UNITS =
   * FT].
   */
  double friction_slope(HydraulicShape *shape, double Q, double depth);
  double n; /**< Manning's n coefficient used in Manning's Equation.*/
};

/**
 * @brief An implementation of the Hazen-Williams Equation to calculate the
 * slope of the Energy Grade Line (EGL) in a Conduit.
 *
 */
class HazenFriction : public FrictionMethod {
public:
  /**
   * @brief Construct a new Hazen Friction object
   *
   * @param C Hazen roughness coefficient.
   */
  HazenFriction(double C);
  /**
   * @brief An implementation of Hazen-Williams Equation to calculate the slope
   * of the Energy Grade Line (EGL) in a Conduit.
   * @param shape The cross-sectional shape of the Conduit.
   * @param Q The flow through the Conduit [UNITS = CFS].
   * @param depth The depth of flow in the Conduit at the calculation point
   * [UNITS = FT].
   * @return double The slope of the EGL at a point along the Conduit[UNITS =
   * FT].
   */
  double friction_slope(HydraulicShape *shape, double Q, double depth);
  double C; /**< C friction coefficient used in the Hazen-Williams Equation.*/
};

/**
 * @brief An implementation of the Darcy-Weisbach Equation to calculate the
 * slope of the Energy Grade Line (EGL) in a Conduit.
 *
 */
class DarcyFriction : public FrictionMethod {
public:
  /**
   * @brief Construct a new Darcy Friction object
   *
   * @param eps Pipe roughness coefficient [UNITS = FT].
   * @param rho Fluid density [UNITS = SLUG/CF].
   * @param mu Fluid dynamic viscosity [UNITS = PSF*S].
   */
  DarcyFriction(double eps, double rho, double mu);
  /**
   * @brief An implementation of Darcy-Weisbach Equation to calculate the slope
   * of the Energy Grade Line (EGL) in a Conduit.
   * @return double The slope of the EGL at a point along the Conduit[UNITS =
   * FT].
   *
   */
  double friction_slope(HydraulicShape *shape, double Q, double depth);
  double eps; /**< epsilon pipe roughness coefficient used in the Darcy-Weisbach
     Equation. Dimension matches that of the Hydraulic Shape [UNITS = FT].*/
  double rho; /**< Density of the fluid [UNITS = SLUG/CF].*/
  double mu;  /**< Dynamic Viscosity of the fluid [UNIT = PSF*S].*/
};

} // namespace hazen

#endif