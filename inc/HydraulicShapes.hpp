#ifndef HYDRAULICSHAPES
#define HYDRAULICSHAPES

// Standard Library
#include <cmath>

namespace hazen {
/**
 * @brief A Hydraulic Shape used by a Hydraulic Link to describe the
 * cross-sectional area of some Hydraulic Component.
 * @details All Hydraulic Shapes provide functions for the top width, wetted
 * perimeter, flow area, hydraulic radius, and hydraulic diameter at a user
 * specified depth.
 *
 */
class HydraulicShape {
public:
  /**
   * @brief The Hydraulic Radius of the Hydraulic Shape.
   * @details Hydraulic Radius = Flow Area / Wetted Perimeter.
   *
   * @param depth The depth of flow to use in the calculation.
   * @return double The Hydraulic Radius at the user specified depth.
   */
  virtual double hydraulic_radius(double depth);
  /**
   * @brief The Hydraulic Diameter of the Hydraulic Shape.
   * @details Hydraulic Diameter = 4 * Hydraulic Radius.
   *
   * @param depth The depth of flow to use in the calculation.
   * @return double The Hydraulic Diameter at the user specified depth.
   */
  virtual double hydraulic_diameter(double depth);
  /**
   * @brief The Top Width of the surface of the water in the cross-section at a
   * user specified depth.
   *
   * @param depth The depth of flow to use in the calculation.
   * @return double The Top Width at the user specified depth. Returns zero if
   * depth is less than or equal to zero. Returns zeros if depth is larger than
   * or equal to the largest vertical dimension of the shape.
   */
  virtual double top_width(double depth) = 0;
  /**
   * @brief The Wetted Perimeter of the cross-section at a user specified depth.
   *
   * @param depth The depth of flow to use in the calculation.
   * @return double The Wetted Perimeter at the user specified depth. Returns
   * zero if depth is less than or equal to zero. Returns the total perimeter of
   * the shape if depth is larger than or equal to the largest vertical
   * dimension of the shape.
   */
  virtual double wetted_perimeter(double depth) = 0;
  /**
   * @brief The Flow Area of the cross-section at a user specified depth.
   *
   * @param depth The depth of flow to use in the calculation.
   * @return double The Flow Area at the user specified depth. Returns
   * zero if depth is less than or equal to zero. Returns the total area of
   * the shape if depth is larger than or equal to the largest vertical
   * dimension of the shape.
   */
  virtual double flow_area(double depth) = 0;
};

/**
 * @brief A Circular cross-section shape.
 *
 */
class Circle : public HydraulicShape {
public:
  Circle(double diameter);
  double top_width(double depth);
  double wetted_perimeter(double depth);
  double flow_area(double depth);
  double height; /**< The diameter of the Circle.*/
};

/**
 * @brief A Rectangular cross-section shape.
 *
 */
class Rectangle : public HydraulicShape {
public:
  Rectangle(double width, double height, bool is_open = true);
  double top_width(double depth);
  double wetted_perimeter(double depth);
  double flow_area(double depth);
  double length; /**< The horizontal length of the Rectangle.*/
  double height; /**< The vertical height of the Rectangle.*/
  bool is_open;  /**< If true, the shape is treated as if it had infinite
                    height.*/
};

/**
 * @brief A V-Notch cross-section shape.
 *
 */
class VNotch : public HydraulicShape {
  VNotch(double top_length, double height, double overflow_length,
         bool is_open = true);
  double top_width(double depth);
  double wetted_perimeter(double depth);
  double flow_area(double depth);
  double top_length; /**< The horizontal length of the top of the v-notch when
                        depth = height.*/
  double overflow_length; /**< The horizontal length to use when depth > height
                             and is_open = true.*/
  double height;          /**< The vertical height of the v-notch.*/
  bool is_open; /**< If true, the v-notch is treated as if it has no top and
                   the overflow length is used for the horizontal dimension with
                   depths larger than the height of the v-notch.*/
};

/**
 * @brief A Trapezoidal cross-section shape.
 *
 */
class Trapezoid : public HydraulicShape {
public:
  Trapezoid(double top_length, double bottom_length, double height,
            double overflow_length, bool is_open = true);
  double top_width(double depth);
  double wetted_perimeter(double depth);
  double flow_area(double depth);
  double bottom_length; /**< The horizontal length of the bottom of the
                           trapezoid. Must be larger than or equal to zero. A
                           zero here represents a V-Notch shape.*/
  double top_length; /**< The horizontal length of the top of the trapezoid when
                        depth = height.*/
  double overflow_length; /**< The horizontal length to use when depth > height
                             and is_open = true.*/
  double height;          /**< The vertical height of the trapezoid.*/
  bool is_open; /**< If true, the trapezoid is treated as if it has no top and
                   the overflow length is used for the horizontal dimension with
                   depths larger than the height of the trapezoid.*/
};

} // namespace hazen

#endif