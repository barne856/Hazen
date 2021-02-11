#ifndef HAZEN_HYDRAULIC_LINKS
#define HAZEN_HYDRAULIC_LINKS

// Standard Library
#include <utility>
#include <vector>

// HAZEN
#include "Hazen/HydraulicNetwork.hpp"

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
              std::pair<Vec<Length>, Vec<Length>> alignment);
  Length velocity_head(Length depth) const;
  Length invert(HydraulicNode *node) const;
  Length length() const;
  Length horizontal_length() const;
  Angle slope(HydraulicNode *node) const;
  Length critical_depth() const;
  Length normal_depth(HydraulicNode *node) const;
  bool is_steep(HydraulicNode *node) const;
  Angle friction_slope(Length depth) const;
  Angle hydrualic_slope_subcritical(Length x, Length h,
                                    HydraulicNode *node) const;
  Angle hydrualic_slope_supercritical(Length x, Length h,
                                      HydraulicNode *node) const;
  /**
   * @brief Calculate the head loss through a conduit or channel of the given
   * cross-sectional shape.
   *
   * @return double The head loss through the conduit or channel from the
   * dowstream node to the upstream node.
   */
  Length head_loss(HydraulicNode *node) override;
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
  Length get_upstream_water_surface(HydraulicNode *node) const;
  Length get_water_surface_subcritical(HydraulicNode *node) const override;
  Length get_water_surface_supercritical(HydraulicNode *node) const;
  std::shared_ptr<HydraulicShape>
      cross_section_shape; /**< The cross-sectional shape of the passage.*/
  std::shared_ptr<FrictionMethod>
      friction_method; /**< The friction method used to calculate the head
                          loss through the passage.*/
  std::pair<Vec<Length>, Vec<Length>>
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
  OpeningLink(std::shared_ptr<HydraulicShape> cross_section_shape,
              Dimensionless Cd, Length elevation,
              Dimensionless percent_open = Dimensionless(1.0));
  /**
   * @brief Calculate the head loss through the opening of the given
   * cross-sectional shape.
   *
   * @return double The head loss through the opening from the
   * dowstream node to the upstream node.
   */
  virtual Length head_loss(HydraulicNode *node) override;
  Length get_water_surface_subcritical(HydraulicNode *node) const override;
  Velocity velocity(Length y, Length d1, Length d2) const;
  quotient_type<Area, Time> flow_step(Length y, Flow Q, Length d1,
                                      Length d2) const;
  Flow flow(Length d1, Length d2) const;

  std::shared_ptr<HydraulicShape>
      cross_section_shape; /**< The cross-sectional shape of the opening.*/
  Dimensionless Cd;        /**< The empirical flow coefficent for the Weir.*/
  Length elevation; /**< The elevation of the lowest point on the weir crest
                       [UNITS = FT].*/
  Dimensionless percent_open; /**< The percent open if this is a gate. 0 to 1*/
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
  MinorLink(Dimensionless K_pos, Dimensionless K_neg,
            std::shared_ptr<HydraulicShape> cross_section_shape, Length invert);
  /**
   * @brief Calculate the head loss due to minor losses in the Hydraulic Link.
   *
   * @return double The minor head loss through the Hydraulic Link.
   * @see HydraulicLink
   * @see HydraulicComponent
   */
  Length head_loss(HydraulicNode *node) override;
  Length get_water_surface_subcritical(HydraulicNode *node) const override;
  Length velocity_head(Length depth) const;

  Dimensionless K_pos; /**< The minor loss coefficent for positive flow */
  Dimensionless K_neg; /**< The minor loss coefficent for negative flow */
  std::shared_ptr<HydraulicShape>
      cross_section_shape; /**< The cross-sectional shape of the opening.*/
  Length invert;
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
  TransitionLink(Dimensionless K,
                 std::pair<std::shared_ptr<HydraulicShape>,
                           std::shared_ptr<HydraulicShape>>
                     shapes,
                 std::pair<Length, Length> inverts);
  /**
   * @brief Calculate the head loss due to the transition of the Hydraulic Link.
   *
   * @return double The transition loss through the Hydraulic Link.
   * @see HydraulicLink
   * @see HydraulicComponent
   */
  Length head_loss(HydraulicNode *node) override;
  Length get_water_surface_subcritical(HydraulicNode *node) const override;
  Length invert(HydraulicNode *node) const;
  HydraulicShape *shape(HydraulicNode *node) const;
  Length velocity_head(HydraulicNode *node, Length depth) const;
  Velocity velocity(HydraulicNode *node, Length depth) const;

  Dimensionless K; /**< The empirical loss coefficent 0 <= K <= 1 */
  std::pair<std::shared_ptr<HydraulicShape>, std::shared_ptr<HydraulicShape>>
      shapes; /**< The cross-sectional shapes of the transition.*/
  std::pair<Length, Length>
      inverts; /**< The invert elevations of the transition shapes.*/
};

class SpecialLossLink : public HydraulicLink {
public:
  Length head_loss(HydraulicNode *node);
  std::vector<Flow> Q_vec;
  std::vector<Length> h_vec;
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