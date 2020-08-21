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
 * @brief Compute the head loss though a passage defined by an alignment with a
 * constant cross-sectional shape.
 * @details A passage may be a closed coundit or an open channel. All bends are
 * assumed to be small and no minor losses are computed as part of the
 * calculation due to bends in the alignment.
 *
 * @param shape The cross-sectional shape of the passage. May be a closed
 * conduit or an open channel.
 * @param friction The friction method used to determine the slope of the energy
 * grade line.
 * @param invert_alignment The alignment of the invert of the passage at the
 * center of the shape starting from the downstream end. Alignments may be
 * vertical.
 * @param Q The flow though the passage.
 * @param E The downstream energy head. This value is used to determine the
 * initial energy head of the passage. If the downstream energy head is less
 * than the energy head due to the brink depth and velocity of the passage, the
 * brink depth energy head is used as the starting energy head for the head loss
 * calculation.
 * @param dd The spatial integration step along the depth of flow perpendicular
 * the the passage alignment.
 * @param ds The spatial integration step along the passage alignment.
 * @param s_start An optional parameter to specify a distance along the
 * alignment to use as the start of the calculation. Default is zero which means
 * the downstream end is used as the starting point for calculations.
 * @param s_final An optional parameter to specify a distance along the
 * alignment to use as the stopping point of the calculation. Default is
 * infinity which means the upstream end is used as the ending point for
 * calculations.
 * @return double The head loss from the given downstream energy head to the
 * upstream energy head of the passage. If the downstream energy head is lower
 * than the starting invert elevation of the passage, this value includes the
 * energy loss due to the downstream energy head being below the starting invert
 * elevation of the alignment.
 */
double passage_loss(HydraulicShape *shape, FrictionMethod *friction,
                    std::vector<vec3> invert_alignment, double Q, double E,
                    double dd, double ds, double s_start = 0.0,
                    double s_final = std::numeric_limits<double>::infinity());
double opening_loss(HydraulicShape *shape, double Cd, double opening_invert,
                    double Q, double E, double dd, double percent_open = 1.0);
double pump_loss(std::vector<std::pair<double, double>> flow_head_curve,
                 double discharge_elevation, double Q, double E);
} // namespace hazen

#endif