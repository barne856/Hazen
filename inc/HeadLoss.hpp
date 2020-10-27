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
 * @brief Compute the total energy head loss through a vertical shaft.
 *
 * @param shape The cross sectional shape of the shaft.
 * @param friction The friction method used.
 * @param up_inv The upstream invert elevation
 * @param dn_inv The downstream invert elevation.
 * @param Q The flow through the shaft.
 * @param h The downstream water surface elevation.
 * @param HGL A reference to the HGL vector where the HGL will be recorded.
 * @return double The total energy head upstream of the shaft.
 */
double vertical_shaft_loss(HydraulicShape *shape, FrictionMethod *friction,
                           double up_inv, double dn_inv, double Q, double h,
                           std::vector<std::pair<double, double>> &HGL);

/**
 * @brief Compute the total energy head upstream of a prismatic channel or
 * conduit using the Gradually Varied Flow (GVF) equation. The function does not
 * use the small angle assumption and so works well with aggressivly sloped
 * pipes. For vertical pipes, use the vertical_shaft_loss() function.
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
 * If a hydraulic jump does not occur, jump_x will be set to NaN.
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
 * @param h The downstream piezometric head. This value is used to determine
 * the downstream control point of the passage. If the downstream piezometric
 * head is less than the elevation of the critical depth, the critical depth is
 * used as the downstream control point for the head loss calculation, otherwise
 * this value is used as the downstream control point [UNITS = FT].
 * @param ds The spatial integration step used in the RK4 method (measured
 * along the pipe length) [UNITS = FT].
 * @param jump_x If a hydraulic jump occurs, the location of the jump as
 * measured horizontally from the downstream end of the pipe will be returned in
 * this variable [UNITS = FT].
 * @param HGL A reference to the HGL vector where the HGL will be recorded.
 * @return double The upstream energy head of the passage or NaN if a
 * hydraulic jump occured [UNITS = FT].
 */
double gvf_backwater_loss(HydraulicShape *shape, FrictionMethod *friction,
                          vec3 up_inv, vec3 dn_inv, double Q, double h,
                          double ds, double &jump_x,
                          std::vector<std::pair<double, double>> &HGL);

/**
 * @brief Compute the total energy head upstream of a prismatic channel or
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
 * @param h The upstream piezometric head. This value is used to determine
 * the upstream control point of the passage. If the upstream piezometric head
 * is greater than the elevation of the critical depth, the critical depth is
 * used as the upstream control point for the head loss calculation, otherwise
 * this value is used as the upstream control point [UNITS = FT].
 * @param ds The spatial integration step used in the RK4 method (measured
 * along the pipe length) [UNITS = FT].
 * @param jump_x The location of a hydraulic jump is there is one. Calculation
 * of the HGL stops at this point.
 * @param HGL A reference to the HGL vector where the HGL will be recorded.
 * @return double The upstream energy head [UNITS = FT].
 */
double gvf_frontwater_loss(HydraulicShape *shape, FrictionMethod *friction,
                           vec3 up_inv, vec3 dn_inv, double Q, double h,
                           double ds, double jump_x,
                           std::vector<std::pair<double, double>> &HGL);

/**
 * @brief Compute the total energy head upstream of an opening with a constant
 * cross-sectional shape and a negligible approach velocity.
 * @details The opening can be a sharp-crested weir/orifice or a broad crested
 * weir/orifice. Both free, submerged, and partially submerged weirs/orfices are
 * handled. Use the approprriate Cd coefficient to describe the expected
 * scenario.
 *
 * The function numerically integrates the flow velocity as determined by
 * Bernoulli's equation as a function of depth from the crest of the opening and
 * multiplies by an empirical flow coefficient to determine the flow. The secant
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
 * @param Q The flow though the opening [UNITS = CFS].
 * @param h The downstream piezometric head. This may be larger, smaller, or
 * equal to the \p opening_invert [UNITS = FT].
 * @param dy The integration step used to numerically intergrate velocity along
 * the cross section of the opening [UNITS = FT].
 * @param HGL A reference to the HGL vector where the HGL will be recorded.
 * @param percent_open The percentage of the opening that is clear of
 * obstruction from a possible gate closed on the opening from above.
 * @return double The total energy head upstream of the opening.
 */
double opening_loss(HydraulicShape *shape, double Cd, double opening_invert,
                    double Q, double h, double dy,
                    std::vector<std::pair<double, double>> &HGL,
                    double percent_open = 1.0);

/**
 * @brief Compute the total energy head upstream of a transition (contraction
 * or expansion) in a conduit or open channel with subcritical flow.
 *
 * @param up_shape The upstream cross sectional shape.
 * @param dn_shape The donwstream cross sectional shape.
 * @param up_inv The upstream invert elevation.
 * @param dn_inv The downstream invert elevation.
 * @param Q The flow through the transition.
 * @param K The empirical loss coefficient 0 <= K <= 1. K ~ 1 for sudden
 * expansion/constraction.
 * @param h The downstream water surface elevation.
 * @param HGL A reference to the HGL vector where the HGL will be recorded.
 * @return double The total energy head upstream of the transition.
 */
double transition_loss(HydraulicShape *up_shape, HydraulicShape *dn_shape,
                       double up_inv, double dn_inv, double Q, double K,
                       double h, std::vector<std::pair<double, double>> &HGL);

/**
 * @brief An enum class to describe the type of benching in a manhole.
 * @details See FHWA HEC-22 Figure 7-6.
 *
 */
enum class BENCH_CONFIGURATION { FLAT, DEPRESSED, HALF, FULL, IMPROVED };

/**
 * @brief Compute the benching coefficient for a manhole.
 *
 * @param D The hydraulic diameter of the outlet pipe.
 * @param d The water depth from the invert of the outlet pipe.
 * @param bench The bench configuration.
 * @return double C_B, the benching coefficient.
 */
double benching_coefficient(BENCH_CONFIGURATION bench, double D, double d);

/**
 * @brief Compute the headloss through a manhole.
 * @details The calculation does not include any entrance or exit losses from
 * pipes. Angled inflow, plunging, and benching are taken into account.
 *
 * @param E_ai Initial estimate of energy level. The total energy head in the
 * manhole after the entrance losses from outgoing pipes (Energy head is
 * approximatley water surface elevation).
 * @param E_i Total energy head downstream of the manhole (use average of
 * downstream values weighted by flow).
 * @param Q_0 Total flow flowing out of the manhole.
 * @param D_0 Outlet hydraulic diameter when flowing full (use maximum of
 * outflowing pipes).
 * @param invert_0 The outlet pipe invert (use minimum of the outlet pipes).
 * @param bench The type of benching in a manhole.
 * @param Q vector of flows from inflow pipes (individual pipe upstream flow)
 * @param invert Invert of the respective upstream pipes
 * @param theta angles, in degrees, between individual upstream pipes and
 * outlet pipe. Elements must have a value between 0 and 180.
 * @param HGL A reference to the HGL vector where the HGL will be recorded.
 * @return double Total energy head in the manhole
 */
double manhole_loss(double E_ai, double E_i, double Q_0, double D_0,
                    double invert_0, BENCH_CONFIGURATION bench,
                    std::vector<double> Q, std::vector<double> invert,
                    std::vector<double> theta,
                    std::vector<std::pair<double, double>> &HGL);

// appoach velocity "u"
// Bar screen velocity "v"
// "n" equal size openings
// shape of one of the "n" openings
// K = 1/0.7, coefficient
// headloss = (v^2-u^2)/(2*0.7*g)

// objective function = f(hl) = (v^2 - u(hl)^2)/(0.7*2*g) - hl = 0
// goal = 0

// upstream channel shape

/**
 * @brief The headloss through a bar screen, assuming the bar screen is clean.
 *
 * @param upstream_channel_shape The shape of the channel upstream of the bar
 * screen.
 * @param channel_invert The upstream channel invert elevation.
 * @param bar_screen_opening_shape The shape of a single opening in the bar
 * screen.
 * @param n The number of openings in the bar screen that flow can travel
 * through.
 * @param Q The flow through the bar screen.
 * @param h The downstream water surface elevation.
 * @return double The energy head upstream of the bar screen.
 */
double bar_screen_loss(HydraulicShape *upstream_channel_shape,
                       double channel_invert,
                       HydraulicShape *bar_screen_opening_shape, int n,
                       double Q, double h);

} // namespace hazen

#endif