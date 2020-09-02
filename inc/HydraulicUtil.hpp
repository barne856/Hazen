#include "FrictionMethods.hpp"
#include "HydraulicShapes.hpp"
#include <functional>
#include <limits>
#include <vector>

namespace hazen {
// Constants -----------------------------
const double g = 32.17405;               // acceleration due to gravity [FT/S^2]
const double rho_water_32 = 62.416 / g;  // density of water at 32F
const double rho_water_40 = 62.432 / g;  // density of water at 40F
const double rho_water_50 = 62.408 / g;  // density of water at 50F
const double rho_water_60 = 62.366 / g;  // density of water at 60F
const double rho_water_70 = 62.300 / g;  // density of water at 70F
const double rho_water_80 = 62.217 / g;  // density of water at 80F
const double rho_water_90 = 62.118 / g;  // density of water at 90F
const double rho_water_100 = 61.998 / g; // density of water at 100F
const double mu_water_32 = 3.732e-5;     // dynamic viscosity of water at 32F.
const double mu_water_40 = 3.228e-5;     // dynamic viscosity of water at 40F.
const double mu_water_50 = 2.730e-5;     // dynamic viscosity of water at 50F.
const double mu_water_60 = 2.344e-5;     // dynamic viscosity of water at 60F.
const double mu_water_70 = 2.034e-5;     // dynamic viscosity of water at 70F.
const double mu_water_80 = 1.791e-5;     // dynamic viscosity of water at 80F.
const double mu_water_90 = 1.500e-5;     // dynamic viscosity of water at 90F.
const double mu_water_100 = 1.423e-5;    // dynamic viscosity of water at 100F.
// ---------------------------------------

/**
 * @brief A simple 3 component structure used for storing vertex information.
 *
 */
struct vec3 {
  double x, y, z;
};

// An n-sized vector of vertices describing an alignment of a passage.
using alignment = std::vector<vec3>;

// Dimensionless numbers -----------------
/**
 * @brief Compute the Reynolds number for the flow.
 *
 * @param rho The fluid density [UNITS = SLUG/CF].
 * @param mu Fluid dynamic viscosity [UNITS = PSF*S].
 * @param V Average fluid velocity [UNITS FPS].
 * @param Dh Passage hydraulic diameter [UNITS = FT].
 * @return double The Reynolds number for the flow.
 */
double reynolds(double rho, double mu, double V, double Dh);

/**
 * @brief Compute the Froude number of the flow.
 * @details
 *
 * @param V Average fluid velocity [UNITS FPS].
 * @param h The hydraulic depth [UNITS = FT].
 * @return double The Froude number of the flow.
 */
double froude(double V, double h);
// ---------------------------------------

/**
 * @brief Compute the darcy friction factor used in the Darcy-Weisbach Equation.
 * @details Valid approximation for all flow regimes. See Bellos, Nalbantis,
 * Tsakiris (2018).
 *
 * @param Re Reynolds Number of the flow
 * @param Dh Hydraulic Diameter of the passage [UNITS = FT].
 * @param eps Passage roughness [UNITS = FT].
 * @return double The darcy friction factor for the passage.
 */
double darcy_friction_factor(double Re, double Dh, double eps);

/**
 * @brief Compute the flow velocity as a function of depth from the surface of
 * an infinite reservoir discharging into the atmosphere.
 *
 * @param h The depth from the surface of the infinite reservoir.
 * @return double The flow velocity at the given depth.
 */
double bernoulli_velocity(double h);

/**
 * @brief Calculate the normal depth (uniform flow depth) at the current flow
 * for the conduit or channel.
 * @details Netwon-Raphson algorithm is used to find the normal depth with
 * tolerance of convergence = 0.0001 * shape height. The initial value of the
 * algorithm is 0.1 * shape height.
 *
 * @param shape The cross-sectional shape of the passage.
 * @param friction The friction method to use in the calculation.
 * @param Q The flow to use in the calculation [UNITS = CFS].
 * @param S The slope of the passage [UNITS = FT/FT].
 *
 * @return double The normal depth at the current flow rate. The absolute
 * value of the signed flow and slope is used in the calculation, NaN is
 * returned if there is no flow, infinity is returned if uniform flow cannot be
 * achieved at the current flow rate.
 */
double normal_depth(HydraulicShape *shape, FrictionMethod *friction, double S,
                    double Q);

/**
 * @brief Compute the critical depth of a passage with constant cross-sectional
 * shape.
 * @details Netwon-Raphson algorithm is used to find the critical depth with
 * tolerance of convergence = 0.0001 * shape height. The initial value of the
 * algorithm is 0.5 * shape height.
 *
 * @param shape The cross-sectional shape of the passage.
 * @param Q The flow to use in the calculation [UNITS = CFS].
 * @return double The critical depth for the passage [UNITS = FT].
 */
double critical_depth(HydraulicShape *shape, double Q);

/**
 * @brief Compute the brink depth of a passage.
 * @details The Brink Depth is the minimum of the uniform (normal) depth and the
 * critical depth of the flow through the passage.
 *
 * @param shape The cross-sectional shape of the passage.
 * @param friction The friction method to use in the calculation.
 * @param S The slope of the passage.
 * @param Q The flow to use in the calculation [UNITS = CFS].
 * @return double The brink depth of the passage [UNITS = FT].
 */
double brink_depth(HydraulicShape *shape, FrictionMethod *friction, double S,
                   double Q);

double clamp(double value, double min, double max);

/**
 * @brief Compute the actual length between two points
 *
 * @param down downstream point
 * @param up upstream point
 * @return double The length between the two points.
 */
double length(vec3 down, vec3 up);

/**
 * @brief Compute the horizontal length between two points.
 *
 * @param down downstream point
 * @param up upstream point
 * @return double The horizontal length between the two points.
 */
double horizontal_length(vec3 down, vec3 up);

/**
 * @brief Compute the slope between and upstream and downstream point.
 * @details The slope is positive when the dowsntream point is lower than the
 * upstream point.
 *
 * @param down The downstream point
 * @param up The upstream point
 * @return double The slope from downstream to upstream.
 */
double slope(vec3 down, vec3 up);

/**
 * @brief Compute the length along a 3-dimensional alignment.
 *
 * @param align The alignment.
 * @return double The true length of the alignment [UNITS = FT].
 */
double alignment_length(alignment &align);

/**
 * @brief Compute the horizontal length of an alignment.
 *
 * @param align The alignment.
 * @return double The horizontal length of the alignment [UNITS = FT].
 */
double alignment_horizontal_length(alignment &align);

/**
 * @brief Find the goal of an objective function using the Secant root
 * finding method.
 *
 * @param goal The goal of the objective function.
 * @param x0 The first initial starting point for the objective function. Should
 * ideally be chosen to lie close to the goal.
 * @param x1 The second initial starting point for the objective function.
 * Should ideally be chosento lie close to the goal.
 * @param TOL The tolerance of the convergence for the solution.
 * @param objective The objective function.
 * @param MAX_ITER The maximum iterations to perform before the function quits
 * trying to refine the estimate.
 * @return double The input to the objective function that results in an
 * approximate value of goal being returned. NaN is returned if the goal cannot
 * be found or if the maximum iterations is exceeded.
 */
double find_goal_secant(double goal, double x0, double x1, double TOL,
                        std::function<double(double)> objective,
                        const int MAX_ITER = 100);

/**
 * @brief Find the goal of an objective function using the Bisection root
 * finding method.
 *
 * @param goal The goal of the objective function.
 * @param a The first initial starting point for the objective function. Should
 * be chosen with "b" to bracket to root. objective(a)*objective(b) must be
 * less than 0.
 * @param b The second initial starting point for the objective function.
 * Should be chosen with "a" to bracket to root. objective(a)*objective(b) must
 * be less than 0.
 * @param TOL The tolerance of the convergence for the solution.
 * @param objective The objective function.
 * @param MAX_ITER The maximum iterations to perform before the function quits
 * trying to refine the estimate.
 * @return double The input to the objective function that results in an
 * approximate value of goal being returned. NaN is returned if the goal cannot
 * be found or if the maximum iterations is exceeded.
 */
double find_goal_bisection(double goal, double a, double b, double TOL,
                           std::function<double(double)> objective,
                           const int MAX_ITER = 100);

/**
 * @brief Compute a single step of the Fourth-order Runge-Kutta numerical
 * integration method.
 *
 * @param F The function F(xi,yi) = dy/dx
 * @param yi The initial value of y used to calculate the next value yi+1.
 * @param xi The value of x used to calculate the next value yi+1.
 * @param dx The integration step used by the RK4 method.
 * @return double yi+1 according to the RK4 method with step \p
 * dx.
 */
double RK4(std::function<double(double, double)> F, double yi, double xi,
           double dx);

/**
 * @brief Linearly interpolate a function of a single variable given as a vector of x,y
 * pairs.
 *
 * @param x The x value to interpolate
 * @return double The interpolated y value at x
 */
double interp_1D(std::vector<std::pair<double, double>>& func, double x);

} // namespace hazen
