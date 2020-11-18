#ifndef HYDRAULICCOMPONENTS
#define HYDRAULICCOMPONENTS

#include "HydraulicLinks.hpp"
#include "HydraulicUtil.hpp"

namespace hazen {
/**
 * @brief A Node Component in a Hydraulic Network.
 *
 */
class Node : public HydraulicComponent {
public:
  Node();
  void add_flow(double Q);
  enum binding_point { NODE = 0 };
};

/**
 * @brief An Outfall Component in a Hydraulic Network.
 *
 */
class Outfall : public HydraulicComponent {
public:
  Outfall(double H);
  double H;
  enum binding_point { NODE = 0 };
};

/**
 * @brief An Opening Component in a Hydraulic Network.
 *
 */
class Opening : public HydraulicComponent {
public:
  Opening(std::shared_ptr<HydraulicShape> opening_shape, double Cd,
          double invert);
  enum binding_point { NODE1 = 0, NODE2 };
};

/**
 * @brief A Passage Component in a Hydraulic Network.
 *
 */
class Passage : public HydraulicComponent {
public:
  Passage(std::shared_ptr<HydraulicShape> cross_section_shape,
          std::shared_ptr<FrictionMethod> friction_method,
          std::vector<vec3> alignment);
  enum binding_point { NODE1 = 0, NODE2 };
};

} // namespace hazen

#endif

// Not yet implemented
// Manhole
// Storage
// Pump
// in-line minor loss (Reducers/meters/Bends/etc.)
// Three-way minor loss (Tee/Wye)
// Four-way minor loss (Cross)
// Baffle
// Granular Filter
// Diffuser
// Slide Gate
// Plug Valve
// Check Valve
// Flume