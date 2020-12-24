#include "Hazen/Hydraulics.hpp"
#include <assert.h>
#include <exception>
#include <iostream>
#include <math.h>

namespace hazen {

// Hydraulic Shape -------------------------------------------------------------
HydraulicShape::HydraulicShape(Length height, bool is_open)
    : height(height.val), is_open(is_open) {
  assert(height.val > 0.0);
}
Length HydraulicShape::hydraulic_radius(Length depth) {
  return depth.val > 0.0 ? flow_area(depth) / wetted_perimeter(depth)
                         : Length(0.0);
}
Length HydraulicShape::hydraulic_diameter(Length depth) {
  return Dimensionless(4.0) * hydraulic_radius(depth);
}
Length HydraulicShape::hydraulic_depth(Length depth) {
  if (!is_open && depth >= height) {
    return Length(std::numeric_limits<double>::infinity());
  }
  return depth.val > 0.0 ? flow_area(depth) / top_width(depth) : Length(0.0);
}
Length HydraulicShape::get_max_depth() {
  return is_open ? Length(std::numeric_limits<double>::infinity()) : height;
}
Length HydraulicShape::get_shape_height() { return height; }
Dimensionless HydraulicShape::froude(Flow Q, Length depth) {
  return depth.val <= 0.0
             ? Dimensionless(std::numeric_limits<double>::infinity())
             : hazen::froude(Q / flow_area(depth), hydraulic_depth(depth));
}
bool HydraulicShape::is_free_surface(Length depth) {
  return !is_open && depth >= height ? false : true;
}

// Circle ----------------------------------------------------------------------
Circle::Circle(Length diameter) : HydraulicShape(diameter) {}
Length Circle::top_width(Length depth) {
  auto root_of = (height * height) / 4.0 - power<2, 1>(depth - height / 2.0);
  return root_of.val <= 0.0 ? Length(0.0) : 2.0 * sqrt(root_of);
}
Length Circle::wetted_perimeter(Length depth) {
  double ratio = 1.0 - 2.0 * depth.val / height.val;
  if (ratio <= -1.0) {
    return Length(M_PI * height.val);
  } else if (ratio >= 1.0) {
    return Length(0.0);
  }
  return Length(acos(ratio) * height.val);
}
Area Circle::flow_area(Length depth) {
  Dimensionless ratio = Dimensionless(1.0) - 2.0 * depth / height;
  if (ratio.val <= -1.0) {
    return Area(M_PI * height * height / 4.0);
  } else if (ratio.val >= 1.0) {
    return Area(0.0);
  }
  double theta = acos(ratio.val);
  double root_of = height.val * depth.val - depth.val * depth.val;
  double root = root_of <= 0.0
                    ? 0.0
                    : (height.val / 2.0 - depth.val) * std::sqrt(root_of);
  return Area(height.val <= 0.0 ? 0.0
                                : theta * height.val * height.val / 4.0 - root);
}

// Rectangle -------------------------------------------------------------------
Rectangle::Rectangle(Length width, Length height, bool is_open)
    : width(width.val), HydraulicShape(height, is_open) {}
Length Rectangle::top_width(Length depth) {
  if (depth.val <= 0.0) {
    return Length(0.0);
  } else if (depth >= height && !is_open) {
    return Length(0.0);
  }
  return width;
}
Length Rectangle::wetted_perimeter(Length depth) {
  if (depth <= Length(0.0)) {
    return Length(0.0);
  } else if (depth >= height && !is_open) {
    return Dimensionless(2.0) * width + Dimensionless(2.0) * height;
  }
  return width + Dimensionless(2.0) * depth;
}
Area Rectangle::flow_area(Length depth) {
  if (depth.val <= 0.0) {
    return Area(0.0);
  } else if (depth >= height && !is_open) {
    return width * height;
  }
  return width * depth;
}

// V-Notch ---------------------------------------------------------------------
VNotch::VNotch(Angle angle, Length overflow_width, Length height, bool is_open)
    : angle(angle.val), overflow_width(overflow_width.val),
      HydraulicShape(height, is_open) {}
Length VNotch::top_width(Length depth) {
  if (depth <= Length(0.0)) {
    return Length(0.0);
  } else if (depth >= height) {
    if (!is_open) {
      return Length(0.0);
    }
    return overflow_width;
  }
  Length top_length = Dimensionless(2.0) * height *
                      Dimensionless(tan(angle.as_radians() / 2.0));
  return top_length * depth / height;
}
Length VNotch::wetted_perimeter(Length depth) {
  if (depth <= Length(0.0)) {
    return Length(0.0);
  } else if (depth >= height) {
    Length top_length = Dimensionless(2.0) * height *
                        Dimensionless(tan(angle.as_radians() / 2.0));
    if (!is_open) {
      return top_length + Dimensionless(2.0) * height /
                              Dimensionless(cos(angle.as_radians() / 2.0));
    }
    return overflow_width - top_length +
           Dimensionless(2.0) * height /
               Dimensionless(cos(angle.as_radians() / 2.0)) +
           Dimensionless(2.0) * (depth - height);
  }
  return Dimensionless(2.0) * depth /
         Dimensionless(cos(angle.as_radians() / 2.0));
}
Area VNotch::flow_area(Length depth) {
  if (depth <= Length(0.0)) {
    return Area(0.0);
  } else if (depth >= height) {
    Length top_length = Dimensionless(2.0) * height *
                        Dimensionless(tan(angle.as_radians() / 2.0));
    if (!is_open) {
      return Dimensionless(0.5) * height * top_length;
    }
    return overflow_width * (depth - height) +
           Dimensionless(0.5) * height * top_length;
  }
  Length top_length =
      Dimensionless(2.0) * depth * Dimensionless(tan(angle.as_radians() / 2.0));
  return Dimensionless(0.5) * depth * top_length;
}

// Trapezoid -------------------------------------------------------------------
Trapezoid::Trapezoid(Length cusp_width, Length bottom_width,
                     Length overflow_width, Length height, bool is_open)
    : cusp_width(cusp_width.val), bottom_width(bottom_width.val),
      overflow_width(overflow_width.val), HydraulicShape(height, is_open) {}
Trapezoid::Trapezoid(Angle side_slope, Length bottom_width, Length height,
                     bool is_open)
    : HydraulicShape(height, is_open), bottom_width(bottom_width.val) {
  cusp_width = bottom_width + 2.0 * side_slope.as_slope() * height;
  overflow_width = cusp_width;
}
Length Trapezoid::top_width(Length depth) {
  if (depth <= Length(0.0)) {
    return Length(0.0);
  } else if (depth >= height) {
    if (!is_open) {
      return Length(0.0);
    }
    return overflow_width;
  }
  return bottom_width + (cusp_width - bottom_width) * (depth / height);
}
Length Trapezoid::wetted_perimeter(Length depth) {
  if (depth <= Length(0.0)) {
    return Length(0.0);
  } else if (depth >= height) {
    Length a = (cusp_width - bottom_width) / Dimensionless(2.0);
    Length side_length = sqrt(a * a + height * height);
    if (!is_open) {
      return bottom_width + cusp_width + Dimensionless(2.0) * side_length;
    }
    return bottom_width + Dimensionless(2.0) * side_length +
           Dimensionless(2.0) * (depth - height) +
           (overflow_width - cusp_width);
  }
  Length t = (depth / height) * (cusp_width - bottom_width) + bottom_width;
  Length a = (t - bottom_width) / Dimensionless(2.0);
  Length side_length = sqrt(a * a + depth * depth);
  return Dimensionless(2.0) * side_length + bottom_width;
}
Area Trapezoid::flow_area(Length depth) {
  if (depth <= Length(0.0)) {
    return Area(0.0);
  } else if (depth >= height) {
    if (!is_open) {
      return height * (cusp_width + bottom_width) / Dimensionless(2.0);
    }
    return (depth - height) * overflow_width +
           height * (cusp_width + bottom_width) / Dimensionless(2.0);
  }
  Length t = (depth / height) * (cusp_width - bottom_width) + bottom_width;
  return depth * (t + bottom_width) / Dimensionless(2.0);
}

} // namespace hazen
