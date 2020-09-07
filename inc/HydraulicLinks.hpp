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
 * @details The cross-sectional shape of the conduit or channel is constant and
 * must be a class derived from HydraulicShape class.
 *
 */
class PassageLink : public HydraulicLink {
public:
  PassageLink();
  /**
   * @brief Calculate the head loss through a conduit or channel of the given
   * cross-sectional shape.
   *
   * @return double The head loss through the conduit or channel from the
   * dowstream node to the upstream node.
   */
  double head_loss() override;
  /**
   * @brief Set the cross sectional shape of the conduit or channel.
   *
   * @param shape The shape to use as the cross-section for the conduit or
   * channel.
   */
  void set_cross_section(std::shared_ptr<HydraulicShape> shape);
  /**
   * @brief Set the friction method of the Conduit used in the head loss
   * calculation.
   *
   * @param friction_method The Friction Method to use.
   */
  void set_friction_method(std::shared_ptr<FrictionMethod> friction_method);
  double get_up_node_velocity();
  double get_dn_node_velocity();
  double get_up_node_area();
  double get_dn_node_area();
  alignment invert_alignment; /**< The physical x,y,z coordinates of the
                                 invert of the passage along the alignment in
                                 space [UNITS = FT,FT,FT].*/
  std::shared_ptr<HydraulicShape>
      cross_section_shape; /**< The cross-sectional shape of the passage.*/
  std::shared_ptr<FrictionMethod>
      friction_method; /**< The friction method used to calculate the head
                          loss through the passage.*/
};

/**
 * @brief The Minor Hydraulic Link will compute the local head loss through the
 * Hydraulic Link due to turbulence and fluid-wall separation in fittings,
 * bends, and other junctions where there is no transition in the passage
 * dimensions.
 * @details Downstream velocity is used in all cases.
 *
 */
class MinorLink : public HydraulicLink {
public:
  MinorLink();
  /**
   * @brief Calculate the head loss due to minor losses in the Hydraulic Link.
   *
   * @return double The minor head loss through the Hydraulic Link.
   * @see HydraulicLink
   * @see HydraulicComponent
   */
  double head_loss() override;
  double K; /**< The minor loss coefficent */
  PassageLink *downstream_passage;
};

/**
 * @brief The Transition Hydraulic Link will compute the head loss due to a
 * transition in passage dimensions.
 *
 * @details Both upstream and downstream velocities are considered in the
 * computation.
 *
 */
class TransitionLink : public HydraulicLink {
public:
  TransitionLink();
  /**
   * @brief Calculate the head loss due to the transition of the Hydraulic Link.
   *
   * @return double The transition loss through the Hydraulic Link.
   * @see HydraulicLink
   * @see HydraulicComponent
   */
  double head_loss() override;
  double K; /**< The local loss coefficent */
  PassageLink *downstream_passage;
  PassageLink *upstream_passage;
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
  OpeningLink();
  /**
   * @brief Calculate the head loss through the opening of the given
   * cross-sectional shape.
   *
   * @return double The head loss through the opening from the
   * dowstream node to the upstream node.
   */
  virtual double head_loss() override;
  /**
   * @brief Set the cross sectional shape of the opening.
   *
   * @param shape The shape to use as the cross-section for the opening.
   */
  void set_cross_section(std::shared_ptr<HydraulicShape> shape);
  double Cd;        /**< The empirical flow coefficent for the Weir.*/
  double elevation; /**< The elevation of the lowest point on the weir crest
                       [UNITS = FT].*/

private:
  std::shared_ptr<HydraulicShape>
      cross_section_shape; /**< The cross-sectional shape of the opening.*/
};

// /**
//  * @brief The Pump Hydraulic Link will compute the head loss through a Pump.
//  * @details The head loss will be negative if the Pump imparts head to the
//  * fluid.
//  *
//  */
// class PumpLink : public HydraulicLink {
// public:
//   PumpLink();
//   /**
//    * @brief Calculate the head loss through the Pump.
//    *
//    * @return double The head loss though the Pump. This will be negative if
//    * the Pump imparts head to the fluid.
//    */
//   double head_loss() override;
//   std::vector<std::pair<double, double>>
//       flow_head_curve; /**< The Pump's flow-head curve. first() is flow in
//                           CFS, second() is head in FT.*/
//   double elevation;    /**< The elevation of the pump discharge.*/
// };
//
// /**
//  * @brief The Slide Gate Hydraulic Link will compute the head loss through a
//  * Slide Gate. If fully closed, this Hydraulic Link is disabled from the
//  * Hydraulic Network.
//  * @details The Slide Gate can have any cross-sectional shape derived from
//  the
//  * Hydraulic Shape class. The Slide Gate can be fully closed, fully open, or
//  * open by a percentage.
//  *
//  */
// class SlideGateLink : public OpeningLink {
// public:
//   SlideGateLink();
//   /**
//    * @brief Calculate the head loss though the Slide Gate.
//    *
//    * @return double The head loss though the Slide Gate from the downstream
//    * Hydraulic Node to the upstream Hydraulic Node. NaN is returned if the
//    * Slide Gate is fully closed.
//    */
//   double head_loss() override;
//   double elevation;    /**< The elevation of the gate invert.*/
//   double percent_open; /**< The percentage that the Slide Gate is open. Can
//                           range from 0 to 1.*/
//   bool is_closed; /**< True if the Slide Gate is fully closed and percent
//   open
//                      is less than or equal to zero.*/
// };
//
// /**
//  * @brief The Stop Valve Hydraulic Link will compute the head loss through a
//  * Stop Valve. If closed, this Hydraulic Link is disabled from the Hydraulic
//  * Network.
//  *
//  */
// class StopValveLink : public LocalLink {
// public:
//   StopValveLink();
//   /**
//    * @brief Calculate the head loss though the Stop Valve Link.
//    *
//    * @return double The head loss through the Stop Valve. NaN is returned
//    * if the Stop Valve is closed.
//    */
//   double head_loss() override;
//   double elevation; /**< The elevation of the stop valve invert.*/
//   bool is_closed;   /**< True is the Stop Valve is closed.*/
// };
//
// /**
//  * @brief The Check Valve Hydraulic Link will compute the head loss though a
//  * Check Valve. Flow can only move from upstream to downstream on this
//  * Hydraulic Link.
//  *
//  */
// class CheckValveLink : public LocalLink {
// public:
//   CheckValveLink();
//   /**
//    * @brief Compute the head loss through th Check Valve.
//    *
//    * @return double The head loss though the Check Valve. NaN is returned if
//    * the flow through the Check Valve is negative.
//    */
//   double head_loss() override;
//   double elevation; /**< The elevation of the check valve invert.*/
// };

// Additional Hydraulic Links not yet implemented:
// 1. Flume Link
// 2. Granular Filter Link
// 3. BarRack/Screen Link
// 4. Diffuser Link

} // namespace hazen

#endif