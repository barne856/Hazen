#ifndef HAZEN_HYDRAULIC_COMPONENTS
#define HAZEN_HYDRAULIC_COMPONENTS

#include "Hazen/HydraulicNetwork.hpp"

namespace hazen {
  using Alignment = typename std::vector<Vec<Length>>;

/**
 * @brief A Node Component in a Hydraulic Network.
 *
 */
class VariableHeadNode : public HydraulicComponent {
public:
  VariableHeadNode();
  void add_flow(Flow Q);
  operator std::shared_ptr<HydraulicNode> &() { return nodes[0]; }
  void bind(std::shared_ptr<HydraulicNode> &node);
  enum binding_point { NODE = 0 };
};

/**
 * @brief An Outfall Component in a Hydraulic Network.
 *
 */
class ConstantHeadNode : public HydraulicComponent {
public:
  ConstantHeadNode(Length H);
  operator std::shared_ptr<HydraulicNode> &() { return nodes[0]; }
  void bind(std::shared_ptr<HydraulicNode> &node);
  enum binding_point { NODE = 0 };
  Length H;
};

/**
 * @brief An Opening Component in a Hydraulic Network.
 *
 */
class Opening : public HydraulicComponent {
public:
  Opening(std::shared_ptr<HydraulicShape> opening_shape, Dimensionless Cd,
          Length invert);
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
          std::vector<Vec<Length>> alignment);
  enum binding_point { NODE1 = 0, NODE2 };
  
};

} // namespace hazen

#endif