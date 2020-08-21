namespace hazen {
/**
 * @brief A simple 3 component structure used for storing vertex information.
 *
 */
struct vec3 {
  double x, y, z;
};
} // namespace hazen

/**
 * @brief Calculate the uniform depth at the current flow for the conduit or
 * channel.
 *
 * @return double The uniform flow depth at the current flow rate. The absolute
 * value of the signed flow is used in the calculation, 0 is returned if there
 * is no flow, infinity is returned if uniform flow cannot be achieved at the
 * current flow rate.
 */
double uniform_depth();

double critical_depth();
double brink_depth();

// Reynolds
// Newtons Method