#ifndef HAZEN_FRICTION_METHODS
#define HAZEN_FRICTION_METHODS

// HAZEN
#include "Hazen/HydraulicShapes.hpp"

namespace hazen {
/**
 * @brief An Abstract Friction Method used to calculate head losses in the
 * Passage Hydraulic Link.
 */
class FrictionMethod {
public:
  /**
   * @brief Compute the slope of the Energy Grade Line (EGL) at one point in
   * a Passage due to irreversible friction losses.
   *
   * @param shape The cross-sectional shape of the Passage.
   * @param Q The flow through the Passage.
   * @param depth The depth of flow in the Passage at the point of interest.
   * @return The slope of the EGL at the point of interest along the
   * Passage.
   */
  virtual Angle friction_slope(HydraulicShape *shape, Flow Q, Length depth) = 0;
};

/**
 * @brief An implementation of Manning's Equation to compute the slope of the
 * Energy Grade Line (EGL) in a Passage.
 *
 */
class ManningsFriction : public FrictionMethod {
public:
  /**
   * @brief Construct a new Mannings Friction object.
   *
   * @param n Manning's n coefficient.
   */
  ManningsFriction(Dimensionless n);
  /**
   * @brief An implementation of Manning's Equation to compute the slope of
   * the Energy Grade Line (EGL) in a Passage.
   * @param shape The cross-sectional shape of the Passage.
   * @param Q The flow through the Passage.
   * @param depth The depth of flow in the Passage at the calculation point.
   * @return The slope of the EGL at a point along the Passage.
   */
  Angle friction_slope(HydraulicShape *shape, Flow Q, Length depth) override;

private:
  Dimensionless n; /**< Manning's n coefficient used in Manning's Equation.*/
};

/**
 * @brief An implementation of the Hazen-Williams Equation to compute the
 * slope of the Energy Grade Line (EGL) in a Passage.
 *
 */
class HazenFriction : public FrictionMethod {
public:
  /**
   * @brief Construct a new Hazen Friction object.
   *
   * @param C Hazen roughness coefficient.
   */
  HazenFriction(Dimensionless C);
  /**
   * @brief An implementation of Hazen-Williams Equation to compute the slope
   * of the Energy Grade Line (EGL) in a Passage.
   * @param shape The cross-sectional shape of the Passage.
   * @param Q The flow through the Passage.
   * @param depth The depth of flow in the Passage at the calculation point.
   * @return double The slope of the EGL at a point along the Passage.
   */
  Angle friction_slope(HydraulicShape *shape, Flow Q, Length depth) override;

private:
  Dimensionless
      C; /**< C friction coefficient used in the Hazen-Williams Equation.*/
};

/**
 * @brief An implementation of the Darcy-Weisbach Equation to compute the
 * slope of the Energy Grade Line (EGL) in a Passage.
 *
 */
class DarcyFriction : public FrictionMethod {
public:
  /**
   * @brief Construct a new Darcy Friction object
   *
   * @param eps Pipe roughness coefficient.
   * @param rho Fluid density.
   * @param mu Fluid dynamic viscosity.
   */
  DarcyFriction(Length eps, Density rho, Dynamic_Viscosity mu);
  /**
   * @brief An implementation of Darcy-Weisbach Equation to compute the slope
   * of the Energy Grade Line (EGL) in a Passage.
   * @return The slope of the EGL at a point along the Passage.
   *
   */
  Angle friction_slope(HydraulicShape *shape, Flow Q, Length depth) override;

private:
  Length
      eps; /**< epsilon, pipe roughness coefficient used in the Darcy-Weisbach
              Equation. Dimension matches that of the Hydraulic Shape.*/
  Density rho;          /**< Density of the fluid.*/
  Dynamic_Viscosity mu; /**< Dynamic Viscosity of the fluid.*/
};

} // namespace hazen

#endif