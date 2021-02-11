#ifndef HAZEN_HYDRAULIC_SHAPES
#define HAZEN_HYDRAULIC_SHAPES

// HAZEN
#include "Hazen/Core.hpp"

// Standard Library
#include <cmath>

namespace hazen {
/**
 * @brief A Hydraulic Shape used by a Hydraulic Link to describe the
 * cross-sectional area of some Hydraulic Component.
 * @details All Hydraulic Shapes provide functions for the top width, wetted
 * perimeter, flow area, hydraulic depth, hydraulic radius, and hydraulic
 * diameter at a user specified flow depth.
 *
 */
class HydraulicShape {
public:
  HydraulicShape(Length height, bool is_open = false);
  /**
   * @brief The Hydraulic Radius of the Hydraulic Shape.
   * @details Hydraulic Radius = Flow Area / Wetted Perimeter.
   *
   * @param depth The depth of flow to use in the calculation.
   * @return The Hydraulic Radius at the user specified depth.
   */
  Length hydraulic_radius(Length depth);
  /**
   * @brief The Hydraulic Diameter of the Hydraulic Shape.
   * @details Hydraulic Diameter = 4 * Hydraulic Radius.
   *
   * @param depth The depth of flow to use in the calculation.
   * @return The Hydraulic Diameter at the user specified depth.
   */
  Length hydraulic_diameter(Length depth);
  /**
   * @brief The Hydraulic Depth of the flow.
   * @details Flow Area / Top Width.
   *
   * @param depth The depth of the flow.
   * @return The Hydraulic Depth of the flow.
   */
  Length hydraulic_depth(Length depth);
  /**
   * @brief Get the max depth of the shape.
   *
   * @return If passage is open, infinity is returned, else height is
   * returned.
   */
  Length get_max_depth();
  /**
   * @brief Get the shape height.
   *
   * @return The shape height.
   */
  Length get_shape_height();
  /**
   * @brief Is the flow free surface or pressure driven?
   *
   * @param depth
   * @return true: The flow is free surface.
   * @return false: The flow is pressure driven.
   */
  bool is_free_surface(Length depth);
  /**
   * @brief The Froude number of the shape with a given flow and depth of flow.
   *
   * @param Q The flow.
   * @param depth The depth of flow.
   * @return The Froude number for the flow.
   */
  Dimensionless froude(Flow Q, Length depth);
  /**
   * @brief The Top Width of the surface of the water in the cross-section at a
   * user specified depth.
   *
   * @param depth The depth of flow to use in the calculation.
   * @return The Top Width at the user specified depth. Returns zero if
   * depth is less than or equal to zero. Returns zeros if depth is larger than
   * or equal to the largest vertical dimension of the shape.
   */
  virtual Length top_width(Length depth) = 0;
  /**
   * @brief The Wetted Perimeter of the cross-section at a user specified depth.
   *
   * @param depth The depth of flow to use in the calculation.
   * @return The Wetted Perimeter at the user specified depth. Returns
   * zero if depth is less than or equal to zero. Returns the total perimeter of
   * the shape if depth is larger than or equal to the largest vertical
   * dimension of the shape.
   */
  virtual Length wetted_perimeter(Length depth) = 0;
  /**
   * @brief The Flow Area of the cross-section at a user specified depth.
   *
   * @param depth The depth of flow to use in the calculation.
   * @return The Flow Area at the user specified depth. Returns
   * zero if depth is less than or equal to zero. Returns the total area of
   * the shape if depth is larger than or equal to the largest vertical
   * dimension of the shape.
   */
  virtual Area flow_area(Length depth) = 0;

protected:
  Length height; /**< The maximum height of the shape when a closed top shape is
                    used.*/
  bool is_open;  /**< If true, the shape is treated as if it had infinite
                      height.*/
};

/**
 * @brief A Circular cross-section shape.
 *
 */
class Circle : public HydraulicShape {
public:
  Circle(Length diameter);
  Length top_width(Length depth) override;
  Length wetted_perimeter(Length depth) override;
  Area flow_area(Length depth) override;
};

/**
 * @brief A Rectangular cross-section shape.
 *
 */
class Rectangle : public HydraulicShape {
public:
  Rectangle(Length width, Length height, bool is_open = true);
  Length top_width(Length depth) override;
  Length wetted_perimeter(Length depth) override;
  Area flow_area(Length depth) override;
  Length get_width();
  Length get_height();

private:
  Length width; /**< The horizontal width of the Rectangle.*/
};

/**
 * @brief A V-Notch cross-section shape.
 *
 */
class VNotch : public HydraulicShape {
  VNotch(Angle angle, Length overflow_width, Length height,
         bool is_open = true);
  Length top_width(Length depth) override;
  Length wetted_perimeter(Length depth) override;
  Area flow_area(Length depth) override;

private:
  Angle angle;           /**< The angle of the V-Notch.*/
  Length overflow_width; /**< The horizontal length to use when depth > height
                             and is_open = true.*/
};

/**
 * @brief A Trapezoidal cross-section shape.
 *
 */
class Trapezoid : public HydraulicShape {
public:
  Trapezoid(Length cusp_width, Length bottom_width, Length overflow_width,
            Length height, bool is_open = true);
  Trapezoid(Angle side_slope, Length bottom_width, Length height,
            bool is_open = true);
  Length top_width(Length depth) override;
  Length wetted_perimeter(Length depth) override;
  Area flow_area(Length depth) override;

private:
  Length cusp_width; /**< The horizontal length of the top of the trapezoid when
                        depth = height.*/
  Length bottom_width; /**< The horizontal length of the bottom of the
                           trapezoid. Must be larger than or equal to zero. A
                           zero here represents a V-Notch shape.*/

  Length overflow_width; /**< The horizontal length to use when depth > height
                             and is_open = true.*/
};

} // namespace hazen

#endif