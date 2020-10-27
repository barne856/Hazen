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
  PassageLink(std::unique_ptr<HydraulicShape> cross_section_shape,
              std::unique_ptr<FrictionMethod> friction_method, vec3 up_inv,
              vec3 dn_inv, double ds);
  /**
   * @brief Calculate the head loss through a conduit or channel of the given
   * cross-sectional shape.
   *
   * @return double The head loss through the conduit or channel from the
   * dowstream node to the upstream node.
   */
  double head_loss() override;
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
  double get_upstream_water_surface();
  std::unique_ptr<HydraulicShape>
      cross_section_shape; /**< The cross-sectional shape of the passage.*/
  std::unique_ptr<FrictionMethod>
      friction_method; /**< The friction method used to calculate the head
                          loss through the passage.*/
  vec3 up_inv;         /**< The upstream invert of the Passage.*/
  vec3 dn_inv;         /**< The downstream invert of the Passage.*/
  double ds;
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
  MinorLink(double K_pos, double K_neg);
  /**
   * @brief Calculate the head loss due to minor losses in the Hydraulic Link.
   *
   * @return double The minor head loss through the Hydraulic Link.
   * @see HydraulicLink
   * @see HydraulicComponent
   */
  double head_loss() override;

private:
  double K_pos; /**< The minor loss coefficent for positive flow */
  double K_neg; /**< The minor loss coefficent for negative flow */
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
  TransitionLink(double K);
  /**
   * @brief Calculate the head loss due to the transition of the Hydraulic Link.
   *
   * @return double The transition loss through the Hydraulic Link.
   * @see HydraulicLink
   * @see HydraulicComponent
   */
  double head_loss() override;

private:
  double K; /**< The empirical loss coefficent 0 <= K <= 1 */
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
  OpeningLink(std::unique_ptr<HydraulicShape> cross_section_shape, double Cd,
              double elevation, double dy, double percent_open = 1.0);
  /**
   * @brief Calculate the head loss through the opening of the given
   * cross-sectional shape.
   *
   * @return double The head loss through the opening from the
   * dowstream node to the upstream node.
   */
  virtual double head_loss() override;

private:
  std::unique_ptr<HydraulicShape>
      cross_section_shape; /**< The cross-sectional shape of the opening.*/
  double Cd;               /**< The empirical flow coefficent for the Weir.*/
  double elevation;    /**< The elevation of the lowest point on the weir crest
                          [UNITS = FT].*/
  double dy;           /**< The integration step along the cross section.*/
  double percent_open; /**< The percent open if this is a gate. 0 to 1*/
};

class ManholeLink : public HydraulicLink {
public:
  ManholeLink(double elevation, BENCH_CONFIGURATION bench_config);
  virtual double head_loss() override;

private:
  double elevation;
  BENCH_CONFIGURATION bench_config;
};

/**
 * @brief A Null Link that always reports zero for it's head loss.
 *
 */
class NullLink : public HydraulicLink {
public:
  virtual double head_loss() override;
};

} // namespace hazen

#endif