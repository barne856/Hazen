#ifndef HAZEN_HYDRAULICS
#define HAZEN_HYDRAULICS
#include "Hazen/FrictionMethods.hpp"
#include "Hazen/HydraulicShapes.hpp"
#include "Hazen/Math.hpp"

#include <limits>
#include <string>
#include <utility>
#include <vector>

namespace hazen {
// Physical Constants ----------------------------------------------------------
// acceleration due to gravity
const Acceleration g = 32.17405_fps2;
// density of water
const std::vector<std::pair<Temperature, Density>> rho_water = {
    {32.0_F, 62.416_pcf / g}, {40.0_F, 62.432_pcf / g},
    {50.0_F, 62.408_pcf / g}, {60.0_F, 62.366_pcf / g},
    {70.0_F, 62.300_pcf / g}, {80.0_F, 62.217_pcf / g},
    {90.0_F, 62.118_pcf / g}, {100.0_F, 61.998_pcf / g}};
// dynamic viscosity of water
const std::vector<std::pair<Temperature, Dynamic_Viscosity>> mu_water = {
    {32.0_F, 3.732e-5_psfs}, {40.0_F, 3.228e-5_psfs}, {50.0_F, 2.730e-5_psfs},
    {60.0_F, 2.344e-5_psfs}, {70.0_F, 2.034e-5_psfs}, {80.0_F, 1.791e-5_psfs},
    {90.0_F, 1.500e-5_psfs}, {100.0_F, 1.423e-5_psfs}};

// Dimensionless numbers -------------------------------------------------------
/**
 * @brief Compute the Reynolds number for the flow.
 *
 * @param rho The fluid density.
 * @param mu Fluid dynamic viscosity.
 * @param V Average fluid velocity.
 * @param Dh Passage hydraulic diameter.
 * @return The Reynolds number for the flow.
 */
Dimensionless reynolds(Density rho, Dynamic_Viscosity mu, Velocity V,
                       Length Dh);

/**
 * @brief Compute the Froude number of the flow.
 *
 * @param V Average fluid velocity.
 * @param h The hydraulic depth.
 * @return The Froude number of the flow.
 */
Dimensionless froude(Velocity V, Length h);

/**
 * @brief Compute the darcy friction factor used in the Darcy-Weisbach Equation
 * for pressure driven flow.
 * @details Valid approximation for all flow regimes. See Bellos, Nalbantis,
 * Tsakiris (2018).
 *
 * @param Re Reynolds Number of the flow.
 * @param Dh Hydraulic Diameter of the passage.
 * @param eps Passage roughness.
 * @return The darcy friction factor for the passage.
 */
Dimensionless darcy_friction_factor_pressure_driven(Dimensionless Re, Length Dh,
                                                    Length eps);

/**
 * @brief Compute the darcy friction factor used in the Darcy-Weisbach Equation
 * for free surface flow.
 * @details Valid approximation for all flow regimes. See Bellos, Nalbantis,
 * Tsakiris (2018).
 *
 * @param Re Reynolds Number of the flow.
 * @param Dh Hydraulic Diameter of the passage.
 * @param eps Passage roughness.
 * @return The darcy friction factor for the passage.
 */
Dimensionless darcy_friction_factor_free_surface(Dimensionless Re, Length Dh,
                                                 Length eps);

// Lambert W Function
/**
 * @brief An approximation of the Lambert W function.
 *
 * @param x the input.
 * @return the output.
 */
double lambert_W(double x);

// Characteristic Lengths ------------------------------------------------------
/**
 * @brief Calculate the normal depth (uniform flow depth) at the current flow
 * for the conduit or channel.
 * @details Netwon-Raphson algorithm is used to find the normal depth with
 * tolerance of convergence = 0.0001 * shape height. The initial value of the
 * algorithm is 0.1 * shape height.
 *
 * @param shape The cross-sectional shape of the passage.
 * @param friction The friction method to use in the calculation.
 * @param Q The flow to use in the calculation.
 * @param S The slope of the passage.
 *
 * @return The normal depth at the current flow rate. The absolute
 * value of the signed flow and slope is used in the calculation, NaN is
 * returned if there is no flow, infinity is returned if uniform flow cannot be
 * achieved at the current flow rate.
 */
Length normal_depth(HydraulicShape *shape, FrictionMethod *friction, Angle S,
                    Flow Q);

/**
 * @brief Compute the critical depth of a passage with constant cross-sectional
 * shape.
 * @details Netwon-Raphson algorithm is used to find the critical depth with
 * tolerance of convergence = 0.0001 * shape height. The initial value of the
 * algorithm is 0.5 * shape height.
 *
 * @param shape The cross-sectional shape of the passage.
 * @param Q The flow to use in the calculation.
 * @return The critical depth for the passage.
 */
Length critical_depth(HydraulicShape *shape, Flow Q);

// Helper Functions ------------------------------------------------------------
template <typename T> using HazenRef = std::shared_ptr<T>;
template <typename T, typename... Args> HazenRef<T> make_ref(Args... args) {
  return std::make_shared<T>(args...);
}
std::vector<Vec<Length>> gen_alignment(Angle slope, Length reach,
                                       Length down_invert);

namespace csv_util {

enum class CSV_UNITS {
  FEET = 0,
  METERS,
};

/**
 * @brief Generate a table for exporting to a csv file from a vector of pairs of
 Lengths.
 *
 * @param values A vector of pairs of doubles.
 * @param label1 The label for the first entry in the pairs.
 * @param label2 The label for the second entry in the pairs.
 * @return std::vector<std::pair<std::string, std::vector<double>>> The
 table
 * for writing to a csv file.
 */
std::vector<std::pair<std::string, std::vector<double>>>
gen_table(std::vector<std::pair<Length, Length>> values, std::string label1,
          std::string label2, CSV_UNITS units);

/**
 * @brief Write a CSV file with double type data in columns and string type
 data
 * in headers.
 *
 * @param filename The output filepath
 * @param dataset The data to write
 */
void write_csv(
    std::string filename,
    std::vector<std::pair<std::string, std::vector<double>>> dataset);
} // namespace csv_util

} // namespace hazen
#endif
