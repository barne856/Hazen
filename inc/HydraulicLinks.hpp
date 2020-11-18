#ifndef HYDRAULICLINKS
#define HYDRAULICLINKS

// Standard Library
#include <utility>
#include <vector>

// HAZEN
#include "FrictionMethods.hpp"
#include "HydraulicNetwork.hpp"
#include "HydraulicShapes.hpp"
#include "HydraulicUtil.hpp"

namespace hazen {

/**
 * @brief The Passage Hydraulic Link will compute the head loss through a closed
 * conduit or an open channel.
 * @details The Passage must be prismatic, having constant cross-sectional area
 * and slope. No local losses are computed, only the friction loss due to the
 * passage walls.
 *
 */
class PassageLink : public HydraulicLink {
public:
  PassageLink(std::shared_ptr<HydraulicShape> cross_section_shape,
              std::shared_ptr<FrictionMethod> friction_method,
              std::pair<vec3, vec3> alignment);
  double velocity_head(double depth);
  double invert(HydraulicNode *node);
  double length();
  double horizontal_length();
  double slope(HydraulicNode *node);
  double critical_depth();
  double normal_depth(HydraulicNode *node);
  bool is_steep(HydraulicNode *node);
  double friction_slope(double depth);
  double hydrualic_slope_subcritical(double x, double h, HydraulicNode *node);
  double hydrualic_slope_supercritical(double x, double h, HydraulicNode *node);
  /**
   * @brief Calculate the head loss through a conduit or channel of the given
   * cross-sectional shape.
   *
   * @return double The head loss through the conduit or channel from the
   * dowstream node to the upstream node.
   */
  double head_loss(HydraulicNode *node) override;
  /**
   * @brief Get the upstream water surface of the passage if the passage is
   * steep.
   * @details The function will traverse the hydraulic network through steep
   * upstream passages recursively until a control point can be established and
   * compute frontwater calculations from that point to determine the upstream
   * water surface elevation of the current passage.
   *
   * @return double The upstream water surface elevation of this passage.
   */
  double get_upstream_water_surface(HydraulicNode *node);
  double get_water_surface_subcritical(HydraulicNode *node) override;
  double get_water_surface_supercritical(HydraulicNode *node);
  std::shared_ptr<HydraulicShape>
      cross_section_shape; /**< The cross-sectional shape of the passage.*/
  std::shared_ptr<FrictionMethod>
      friction_method; /**< The friction method used to calculate the head
                          loss through the passage.*/
  std::pair<vec3, vec3>
      alignment; /**< The alignment of the invert of the Passage.*/
  bool is_supercritical = false;
};

/**
 * @brief The Opening Hydraulic Link will compute the head loss through an
 * Orifice or Weir.
 * @details Both sharp crested and broad crested weirs/orifices are
 * allowed. Both submerged and unsubmerged weirs/orifices are allowed.
 *
 */
class OpeningLink : public HydraulicLink {
public:
  OpeningLink(std::shared_ptr<HydraulicShape> cross_section_shape, double Cd,
              double elevation, double percent_open = 1.0);
  /**
   * @brief Calculate the head loss through the opening of the given
   * cross-sectional shape.
   *
   * @return double The head loss through the opening from the
   * dowstream node to the upstream node.
   */
  virtual double head_loss(HydraulicNode *node) override;
  double get_water_surface_subcritical(HydraulicNode *node) override;
  double velocity(double y, double d1, double d2);
  double flow_step(double y, double Q, double d1, double d2);
  double flow(double d1, double d2);

  std::shared_ptr<HydraulicShape>
      cross_section_shape; /**< The cross-sectional shape of the opening.*/
  double Cd;               /**< The empirical flow coefficent for the Weir.*/
  double elevation;    /**< The elevation of the lowest point on the weir crest
                          [UNITS = FT].*/
  double percent_open; /**< The percent open if this is a gate. 0 to 1*/
};

/**
 * @brief The Minor Hydraulic Link will compute the local head loss through the
 * Hydraulic Link due to turbulence and fluid-wall separation in fittings,
 * bends, and other junctions where there is no transition in the passage
 * dimensions.
 * @details The Minor Link must be between exactally two Passage Links.
 * @details Downstream velocity is used in all cases.
 *
 */
class MinorLink : public HydraulicLink {
public:
  MinorLink(double K_pos, double K_neg,
            std::shared_ptr<HydraulicShape> cross_section_shape, double invert);
  /**
   * @brief Calculate the head loss due to minor losses in the Hydraulic Link.
   *
   * @return double The minor head loss through the Hydraulic Link.
   * @see HydraulicLink
   * @see HydraulicComponent
   */
  double head_loss(HydraulicNode *node) override;
  double get_water_surface_subcritical(HydraulicNode *node) override;
  double velocity_head(double depth);

  double K_pos; /**< The minor loss coefficent for positive flow */
  double K_neg; /**< The minor loss coefficent for negative flow */
  std::shared_ptr<HydraulicShape>
      cross_section_shape; /**< The cross-sectional shape of the opening.*/
  double invert;
};

/**
 * @brief The Transition Hydraulic Link will compute the head loss due to a
 * transition in passage dimensions using the Borda-Carnot equation.
 * @details
 *
 * @details Both upstream and downstream velocities are considered in the
 * computation.
 *
 */
class TransitionLink : public HydraulicLink {
public:
  TransitionLink(double K,
                 std::pair<std::shared_ptr<HydraulicShape>,
                           std::shared_ptr<HydraulicShape>>
                     shapes,
                 std::pair<double, double> inverts);
  /**
   * @brief Calculate the head loss due to the transition of the Hydraulic Link.
   *
   * @return double The transition loss through the Hydraulic Link.
   * @see HydraulicLink
   * @see HydraulicComponent
   */
  double head_loss(HydraulicNode *node) override;
  double get_water_surface_subcritical(HydraulicNode *node) override;
  double invert(HydraulicNode *node);
  HydraulicShape *shape(HydraulicNode *node);
  double velocity_head(HydraulicNode *node, double depth);
  double velocity(HydraulicNode *node, double depth);

  double K; /**< The empirical loss coefficent 0 <= K <= 1 */
  std::pair<std::shared_ptr<HydraulicShape>, std::shared_ptr<HydraulicShape>>
      shapes; /**< The cross-sectional shapes of the transition.*/
  std::pair<double, double>
      inverts; /**< The invert elevations of the transition shapes.*/
};

// class ManholeLink : public HydraulicLink {
// public:
//   ManholeLink(double elevation, BENCH_CONFIGURATION bench_config);
//   virtual double head_loss(HydraulicNode *node) override;
//
// private:
//   double elevation;
//   BENCH_CONFIGURATION bench_config;
// };
//
// /**
//  * @brief A Null Link that always reports zero for it's head loss.
//  *
//  */
// class NullLink : public HydraulicLink {
// public:
//   virtual double head_loss(HydraulicNode *node, double h) override;
// };

} // namespace hazen

#endif