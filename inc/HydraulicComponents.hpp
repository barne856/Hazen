#ifndef HYDRAULICCOMPONENTS
#define HYDRAULICCOMPONENTS

#include "HeadLoss.hpp"
#include "HydraulicLinks.hpp"
#include "HydraulicNetwork.hpp"

namespace hazen {
/**
 * @brief An Outfall Component in a Hydraulic Network
 * 
 */
class Outfall : public HydraulicComponent {
public:
  Outfall();
  virtual void bind(std::pair<HydraulicComponent *, unsigned int> bind_point,
                    unsigned int binding_index);

protected:
  virtual HydraulicNode *get_binding_node(unsigned int binding_index);

private:
  std::shared_ptr<HydraulicNode> outfall_node;
};

/**
 * @brief Manhole Component in a Hydraulic Network
 *
 */
class Manhole : public HydraulicComponent {
public:
  Manhole(double invert, BENCH_CONFIGURATION bench_config);
  virtual void bind(std::pair<HydraulicComponent *, unsigned int> bind_point,
                    unsigned int binding_index);

protected:
  virtual HydraulicNode *get_binding_node(unsigned int binding_index);

private:
  double invert;
  BENCH_CONFIGURATION bench_config;
  std::shared_ptr<ManholeLink> manhole_link;
  std::shared_ptr<HydraulicNode> up_node;
  std::shared_ptr<HydraulicNode> dn_node;
};

class Passage : public HydraulicComponent {
public:
  Passage();
  virtual void bind(std::pair<HydraulicComponent *, unsigned int> bind_point,
                    unsigned int binding_index);

protected:
  virtual HydraulicNode *get_binding_node(unsigned int binding_index);

private:
};

/**
 * @brief Storage Component in a Hydraulic Netowork
 * @details When a Passage Link connects to this component, a transition link is
 * placed between the passage and the interior node with loss coefficient 1.
 * Otherwise, a Null Link is used to connect the components.
 *
 */
class Storage : public HydraulicComponent {
public:
  Storage(std::vector<std::pair<double, double>> storage_curve,
          double elevation);
  virtual void bind(std::pair<HydraulicComponent *, unsigned int> bind_point,
                    unsigned int binding_index);

protected:
  virtual HydraulicNode *get_binding_node(unsigned int binding_index);

private:
  std::vector<std::pair<double, double>>
      storage_curve; /**< The storage curve for this type of storage, a vector
                        of pairs (depth, surface area)*/
  double elevation;  /**< The elevation of the bottom of the storage area.*/
  std::shared_ptr<HydraulicNode> storage_node;
};

// Opening


} // namespace hazen

#endif

// Not yet implemented
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