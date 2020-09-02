#ifndef HEADLOSS
#define HEADLOSS

// Standard Library
#include <limits>
#include <utility>
#include <vector>

// HAZEN
#include "FrictionMethods.hpp"
#include "HydraulicShapes.hpp"
#include "HydraulicUtil.hpp"

// Head Loss APIs
namespace hazen {
/**
 * @brief Compute the head loss though a fitting or other junction due to
 * turbulence and fluid-wall separation.
 *
 * @param V The velocity used in the head loss calculation [UNITS = FPS].
 * @param K The minor loss coefficient used in the head loss calculation.
 * @return double The minor head loss [UNITS = FT].
 */
double minor_loss(double V, double K);

/**
 * @brief Compute the piezometric head loss through a prismatic channel or
 * conduit using the Gradually Varied Flow (GVF) equation. The function does not
 * use the small angle assumption and so works well with aggressivly sloped
 * pipes and even vertical pipes.
 * @details This function integrates the GVF equation starting from a downstream
 * control point using the Fourth Order Runge-Kutta (RK4) numerical integration
 * method. This function should only be used with non-steep passages where the
 * normal depth is less than the critical depth or with steep passages where the
 * downstream end is submerged above the critical depth.
 *
 * If this function is used with a steep passage where the downstream end is
 * submerged above the critical depth and a hydraulic jump occurs, the function
 * will return NaN and set the variable \p jump_x to the horizontal distance of
 * the jump as measured from the downstream end of the passage.
 *
 * @param shape The cross-sectional shape of the passage. May be a closed
 * conduit or an open channel.
 * @param friction The friction method used to determine the slope of the
 * energy grade line.
 * @param up_inv The upstream invert coordinates of the passage [UNITS =
 * FT,FT,FT].
 * @param dn_inv The downstream invert coordinates of the passage [UNITS =
 * FT,FT,FT].
 * @param Q The flow though the passage must be positive [UNITS = CFS].
 * @param H The downstream piezometric head. This value is used to determine
 * the downstream control point of the passage. If the downstream piezometric
 * head is less than the elevation of the critical depth, the critical depth is
 * used as the downstream control point for the head loss calculation, otherwise
 * this value is used as the downstream control point [UNITS = FT].
 * @param dx The spatial integration step used in the RK4 method (measured
 * horizontally) [UNITS = FT].
 * @param jump_x If a hydraulic jump occurs, the location of the jump as
 * measured horizontally from the downstream end of the pipe will be returned in
 * this variable [UNITS = FT].
 * @param x_start An optional parameter to specify a horizontal distance along
 * the alignment to use as the start of the calculation. Default is zero which
 * means the downstream end is used as the starting point for calculations
 * [UNITS = FT].
 * @param x_final An optional parameter to specify a horizontal distance along
 * the alignment to use as the stopping point of the calculation. Default is
 * infinity which means the upstream end is used as the ending point for
 * calculations [UNITS = FT].
 * @return double The upstream piezometric head of the passage or NaN if a
 * hydraulic jump occured [UNITS = FT].
 */
double
gvf_backwater_loss(HydraulicShape *shape, FrictionMethod *friction, vec3 up_inv,
                   vec3 dn_inv, double Q, double H, double dx, double &jump_x,
                   double x_start = 0.0,
                   double x_final = std::numeric_limits<double>::infinity());

/**
 * @brief Compute the piezometric head loss through a prismatic channel or
 * conduit using the Gradually Varied Flow (GVF) equation.
 * @details This function integrates the GVF equation starting from an upstream
 * control point using the Fourth Order Runge-Kutta (RK4) numerical integration
 * method. This function should only be used with steep passages where the
 * normal depth is less than the critical depth.
 *
 * @param shape The cross-sectional shape of the passage. May be a closed
 * conduit or an open channel.
 * @param friction The friction method used to determine the slope of the
 * energy grade line.
 * @param up_inv The upstream invert coordinates of the passage [UNITS =
 * FT,FT,FT].
 * @param dn_inv The downstream invert coordinates of the passage [UNITS =
 * FT,FT,FT].
 * @param Q The flow though the passage must be positive [UNITS = CFS].
 * @param H The upstream piezometric head. This value is used to determine
 * the upstream control point of the passage. If the upstream piezometric head
 * is greater than the elevation of the critical depth, the critical depth is
 * used as the upstream control point for the head loss calculation, otherwise
 * this value is used as the upstream control point [UNITS = FT].
 * @param dx The spatial integration step used in the RK4 method (measured
 * horizontally) [UNITS = FT].
 * @param x_start An optional parameter to specify a horizontal distance along
 * the alignment to use as the start of the calculation. Default is zero which
 * means the upstream end is used as the starting point for calculations
 * [UNITS = FT].
 * @param x_final An optional parameter to specify a horizontal distance along
 * the alignment to use as the stopping point of the calculation. Default is
 * infinity which means the downstream end is used as the ending point for
 * calculations [UNITS = FT].
 * @return double The upstream piezometric head of the passage or NaN if a
 * hydraulic jump occured [UNITS = FT].
 */
double
gvf_frontwater_loss(HydraulicShape *shape, FrictionMethod *friction,
                    vec3 up_inv, vec3 dn_inv, double Q, double H, double dx,
                    double x_start = 0.0,
                    double x_final = std::numeric_limits<double>::infinity());

/**
 * @brief Compute the head loss through an opening with a constant
 * cross-sectional shape and an approach velocity.
 * @details The opening can be a sharp-crested weir/orifice or a broad crested
 * weir/orifice. Both free, submerged, and partially submerged weirs/orfices are
 * handled.
 *
 * The function numerically integrates the flow velocity as determined by
 * Bernoulli's equation as a function of depth from the crest of the opening and
 * multiplies by an empirical flow coefficent to determine the flow. Newton's
 * method is used to solve for the head loss needed for the flow though the
 * opening to equal \p Q.
 *
 * @param shape The cross-sectional shape of the opening.
 * @param Cd The weir/orifice discharge coefficient, empirically determined flow
 * coefficient depending on the type of weir/orifice and height of the crest
 * above the channel bottom. Typical values for free flowing sharp crested weir
 * with a vented nappe and a large height above the channel bottom are below:
 * Shape         | Cd
 * ------------- | -------------
 * Rectangular   | 0.62
 * 90 V-Notch    | 0.58
 * @param opening_invert The invert elevation of the opening [UNITS = FT].
 * @param V The approach velocity. This is the average approach velocity of the
 * fluid upstream of the opening [UNITS = FPS].
 * @param Q The flow though the opening [UNITS = CFS].
 * @param H The downstream piezometric head. This may be larger, smaller, or
 * equal to the \p opening_invert [UNITS = FT].
 * @param dy The integration step used to numerically intergrate velocity along
 * the cross section of the opening [UNITS = FT].
 * @param percent_open The percentage of the opening that is clear of
 * obstruction from a possible gate closed on the opening from above.
 * @return double The head loss from the given downstream energy head to the
 * upstream energy head of the opening. If the downstream energy head is lower
 * than the invert elevation of the opening, this value includes the
 * energy loss due to the downstream energy head being below the opening.
 */
double opening_loss(HydraulicShape *shape, double Cd, double opening_invert,
                    double Q, double H, double dy, double percent_open = 1.0);
/**
 * @brief Compute the head loss through a Pump.
 * @details The head loss will be negative if the Pump imparts energy to the
 * fluid.
 *
 * @param flow_head_curve The Pump's flow-head curve. first() is flow in CFS,
 * second() is head in FT.
 * @param discharge_elevation The Pump's discharge elevation. This is the
 * discharge of the actual Pump, not any piping that may be connected to the
 * Pump discharge [UNITS = FT].
 * @param Q The flow through the Pump [UNITS = CFS].
 * @param E The downstream energy head of the Pump [UNITS = FT].
 * @return double The head loss due to the Pump.
 */
double pump_loss(std::vector<std::pair<double, double>> flow_head_curve,
                 double discharge_elevation, double Q, double E);
} // namespace hazen

#endif