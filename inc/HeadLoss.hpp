#ifndef HEADLOSS
#define HEADLOSS

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

// forward declaration
struct vec3;

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
double vertical_drop_loss(PassageLink *link, HydraulicNode *node);

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
double vertical_rise_loss(PassageLink *link, HydraulicNode *node);

double backwater_loss(PassageLink *link, HydraulicNode *node);
double frontwater_loss(PassageLink *link, HydraulicNode *node,
                       double jump_x = 0.0);
double opening_loss(OpeningLink *link, HydraulicNode *node);
double transition_loss(TransitionLink *link, HydraulicNode *node);

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
 * @param D_0 height of outlet opening (use maximum of
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

/**
 * @brief Compute the energy head upstream of a bar screen, assuming the bar
 * screen is clean.
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