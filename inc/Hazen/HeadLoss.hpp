#ifndef HAZEN_HEAD_LOSS
#define HAZEN_HEAD_LOSS

#include "Hazen/HydraulicLinks.hpp"

// Standard Library
#include <limits>
#include <utility>
#include <vector>

// Head Loss APIs
namespace hazen {

// HAZEN
class PassageLink;
class HydraulicNode;
class OpeningLink;
class TransitionLink;
class HydraulicShape;

/**
 * @brief Compute the total energy head at the upstream node of a PassageLink
 * that describes a vertical drop.
 *
 * @details The flow in the PassageLink is assumed to be directed towards the
 * given \p node. The elevation of the invert at this \p node is assumed to be
 * lower than the elevation of the invert at the upstream node. The (x,y)
 * coordinates of the PassageLink's inverts are assumed to be equal, thus the
 * passage should be vertical and discharging downward.
 *
 * @param link The PassageLink
 * @param node The downstream HydraulicNode
 * @return double The energy head at the upstream Hydraulic Node opposite of \p
 * node.
 */
Length vertical_drop_loss(PassageLink *link, HydraulicNode *node);

/**
 * @brief Compute the total energy head at the upstream node of a PassageLink
 * that describes a vertical rise.
 *
 * @details The flow in the PassageLink is assumed to be directed towards the
 * given \p node. The elevation of the invert at this \p node is assumed to be
 * higher than the elevation of the invert at the upstream node. The (x,y)
 * coordinates of the PassageLink's inverts are assumed to be equal, thus the
 * passage should be vertical and discharging upward.
 *
 * @param link The PassageLink
 * @param node The downstream HydraulicNode
 * @return double The energy head at the upstream Hydraulic Node opposite of \p
 * node.
 */
Length vertical_rise_loss(PassageLink *link, HydraulicNode *node);

/**
 * @brief Compute the total energy head at the upstream node of a PassageLink
 * that is backwatered or has a mild slope.
 * @details This computation solves the gradually varied flow (GVF) equations
 * using an adaptive step size 5th order accurate Runge-Kutta method. Only
 * subcritical flow is computed, thus if a hydraulic jump occurs or the head
 * loss is determined by an upstream control point, the function returns NaN and
 * records the HGL profile up the the point of supercritical flow.
 *
 * @param link The PassageLink
 * @param node The downstream HydraulicNode
 * @return double The total energy head at the upstream HydraulicNode opposite
 * of \p node.
 */

Length backwater_loss(PassageLink *link, HydraulicNode *node);
/**
 * @brief Compute the total energy head at the upstream node of a PassageLink
 * that has a steep slope and has its head loss determined by an upstream
 * control point.
 * @details This computation solves the gradually varied flow (GVF) equations
 * using an adaptive step size 5th order accurate Runge-Kutta method. Only
 * supercritical flow is computed, thus if a hydraulic jump occured or the head
 * downstream, the calculation will stop at that point ( See \p jump_x ).
 *
 * @param link The PassageLink
 * @param node The downstream HydraulicNode
 * @param jump_x The distance along the horizontal direction of the hydraulic
 * jump if it was previously copmuted using the backwater_loss() function.
 * @return double The total energy head at the upstream HydraulicNode opposite
 * of \p node.
 */
Length frontwater_loss(PassageLink *link, HydraulicNode *node,
                       Length jump_x = Length(0.0));

/**
 * @brief Compute the total energy head at the upstream node of an OpeningLink.
 * @details An OpeningLink can be used to model an orifice or a weir. It is
 * assumed that there will be no approach velocity and the upstream total energy
 * head is the same as the upstream water surface elevation.
 *
 * @param link The OpeningLink
 * @param node The downstream HydraulicNode
 * @return double The total energy head at the upstream HydraulicNode opposite
 * of \p node.
 */
Length opening_loss(OpeningLink *link, HydraulicNode *node);

/**
 * @brief Compute the total energy head at the upstream node of a
 * TransitionLink.
 * @details The Borda-Carnot Equation is used to model gradual and aprupt
 * transitions (expansion and contraction) in open channel and closed conduit
 * flow.
 *
 * @param link The TransitionLink.
 * @param node The upstream HydraulicNode.
 * @return double The total energy head at the upstream HydraulicNode opposite
 * of \p node.
 */
Length transition_loss(TransitionLink *link, HydraulicNode *node);

// /**
//  * @brief An enum class to describe the type of benching in a manhole.
//  * @details See FHWA HEC-22 Figure 7-6.
//  *
//  */
// enum class BENCH_CONFIGURATION { FLAT, DEPRESSED, HALF, FULL, IMPROVED };
//
// /**
//  * @brief Compute the benching coefficient for a manhole.
//  *
//  * @param D The hydraulic diameter of the outlet pipe.
//  * @param d The water depth from the invert of the outlet pipe.
//  * @param bench The bench configuration.
//  * @return double C_B, the benching coefficient.
//  */
// double benching_coefficient(BENCH_CONFIGURATION bench, double D, double d);
//
// /**
//  * @brief Compute the headloss through a manhole.
//  * @details The calculation does not include any entrance or exit losses from
//  * pipes. Angled inflow, plunging, and benching are taken into account.
//  *
//  * @param E_ai Initial estimate of energy level. The total energy head in the
//  * manhole after the entrance losses from outgoing pipes (Energy head is
//  * approximatley water surface elevation).
//  * @param E_i Total energy head downstream of the manhole (use average of
//  * downstream values weighted by flow).
//  * @param Q_0 Total flow flowing out of the manhole.
//  * @param D_0 height of outlet opening (use maximum of
//  * outflowing pipes).
//  * @param invert_0 The outlet pipe invert (use minimum of the outlet pipes).
//  * @param bench The type of benching in a manhole.
//  * @param Q vector of flows from inflow pipes (individual pipe upstream flow)
//  * @param invert Invert of the respective upstream pipes
//  * @param theta angles, in degrees, between individual upstream pipes and
//  * outlet pipe. Elements must have a value between 0 and 180.
//  * @param HGL A reference to the HGL vector where the HGL will be recorded.
//  * @return double Total energy head in the manhole
//  */
// double manhole_loss(double E_ai, double E_i, double Q_0, double D_0,
//                     double invert_0, BENCH_CONFIGURATION bench,
//                     std::vector<double> Q, std::vector<double> invert,
//                     std::vector<double> theta,
//                     std::vector<std::pair<double, double>> &HGL);
//
// /**
//  * @brief Compute the energy head upstream of a bar screen, assuming the bar
//  * screen is clean.
//  *
//  * @param upstream_channel_shape The shape of the channel upstream of the bar
//  * screen.
//  * @param channel_invert The upstream channel invert elevation.
//  * @param bar_screen_opening_shape The shape of a single opening in the bar
//  * screen.
//  * @param n The number of openings in the bar screen that flow can travel
//  * through.
//  * @param Q The flow through the bar screen.
//  * @param h The downstream water surface elevation.
//  * @return double The energy head upstream of the bar screen.
//  */
// double bar_screen_loss(HydraulicShape *upstream_channel_shape,
//                        double channel_invert,
//                        HydraulicShape *bar_screen_opening_shape, int n,
//                        double Q, double h);
//
} // namespace hazen

#endif