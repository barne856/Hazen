#include "HydraulicShapes.hpp"
#include "HydraulicUtil.hpp"
#include <exception>
#include <iostream>
#define _USE_MATH_DEFINES
#include <assert.h>
#include <math.h>

namespace hazen {

// Hydraulic Shape -------------------------------------------------------------
HydraulicShape::HydraulicShape(double height, bool is_open)
    : height(height), is_open(is_open) {
  assert(height > 0.0);
}
double HydraulicShape::hydraulic_radius(double depth) {
  return depth > 0.0 ? flow_area(depth) / wetted_perimeter(depth) : 0.0;
}
double HydraulicShape::hydraulic_diameter(double depth) {
  return 4.0 * hydraulic_radius(depth);
}
double HydraulicShape::hydraulic_depth(double depth) {
  if (!is_open && depth >= height) {
    return std::numeric_limits<double>::infinity();
  }
  return depth > 0.0 ? flow_area(depth) / top_width(depth) : 0.0;
}
double HydraulicShape::get_max_depth() {
  return is_open ? std::numeric_limits<double>::infinity() : height;
}
double HydraulicShape::get_shape_height() { return height; }
double HydraulicShape::froude(double Q, double depth) {
  return depth <= 0.0
             ? std::numeric_limits<double>::infinity()
             : hazen::froude(Q / flow_area(depth), hydraulic_depth(depth));
}
bool HydraulicShape::is_free_surface(double depth) {
  return !is_open && depth >= height ? false : true;
}

// Circle ----------------------------------------------------------------------
Circle::Circle(double diameter) : HydraulicShape(diameter) {}
double Circle::top_width(double depth) {
  double root_of = (height * height) / 4.0 - pow((depth - height / 2.0), 2.0);
  return root_of <= 0.0 ? 0.0 : 2.0 * sqrt(root_of);
}
double Circle::wetted_perimeter(double depth) {
  double ratio = 1.0 - 2.0 * depth / height;
  if (ratio <= -1.0) {
    return M_PI * height;
  } else if (ratio >= 1.0) {
    return 0.0;
  }
  return acos(ratio) * height;
}
double Circle::flow_area(double depth) {
  double ratio = 1.0 - 2.0 * depth / height;
  if (ratio <= -1.0) {
    return M_PI * height * height / 4.0;
  } else if (ratio >= 1.0) {
    return 0.0;
  }
  double theta = acos(ratio);
  double root_of = height * depth - depth * depth;
  double root = root_of <= 0.0 ? 0.0 : (height / 2.0 - depth) * sqrt(root_of);
  return height <= 0.0 ? 0.0 : theta * height * height / 4.0 - root;
}

// Rectangle -------------------------------------------------------------------
Rectangle::Rectangle(double width, double height, bool is_open)
    : width(width), HydraulicShape(height, is_open) {}
double Rectangle::top_width(double depth) {
  if (depth <= 0.0) {
    return 0.0;
  } else if (depth >= height && !is_open) {
    return 0.0;
  }
  return width;
}
double Rectangle::wetted_perimeter(double depth) {
  if (depth <= 0.0) {
    return 0.0;
  } else if (depth >= height && !is_open) {
    return 2.0 * width + 2.0 * height;
  }
  return width + 2.0 * depth;
}
double Rectangle::flow_area(double depth) {
  if (depth <= 0.0) {
    return 0.0;
  } else if (depth >= height && !is_open) {
    return width * height;
  }
  return width * depth;
}

// V-Notch ---------------------------------------------------------------------
VNotch::VNotch(double angle, double overflow_width, double height, bool is_open)
    : angle(angle), overflow_width(overflow_width),
      HydraulicShape(height, is_open) {}
double VNotch::top_width(double depth) {
  if (depth <= 0.0) {
    return 0.0;
  } else if (depth >= height) {
    if (!is_open) {
      return 0.0;
    }
    return overflow_width;
  }
  double rads = angle * M_PI / 180.0;
  double top_length = 2.0 * height * tan(rads / 2.0);
  return top_length * depth / height;
}
double VNotch::wetted_perimeter(double depth) {
  if (depth <= 0.0) {
    return 0.0;
  } else if (depth >= height) {
    double rads = angle * M_PI / 180.0;
    double top_length = 2.0 * height * tan(rads / 2.0);
    if (!is_open) {
      return top_length + 2.0 * height / cos(rads / 2.0);
    }
    return overflow_width - top_length + 2.0 * height / cos(rads / 2.0) +
           2.0 * (depth - height);
  }
  double rads = angle * M_PI / 180.0;
  return 2.0 * depth / cos(rads / 2.0);
}
double VNotch::flow_area(double depth) {
  if (depth <= 0.0) {
    return 0.0;
  } else if (depth >= height) {
    double rads = angle * M_PI / 180.0;
    double top_length = 2.0 * height * tan(rads / 2.0);
    if (!is_open) {
      return 0.5 * height * top_length;
    }
    return overflow_width * (depth - height) + 0.5 * height * top_length;
  }
  double rads = angle * M_PI / 180.0;
  double top_length = 2.0 * depth * tan(rads / 2.0);
  return 0.5 * depth * top_length;
}

// Trapezoid -------------------------------------------------------------------
Cipolletti::Cipolletti(double top_width, double bottom_width,
                       double overflow_width, double height, bool is_open)
    : cusp_width(cusp_width), bottom_width(bottom_width),
      overflow_width(overflow_width), HydraulicShape(height, is_open) {}
double Cipolletti::top_width(double depth) {
  if (depth <= 0.0) {
    return 0.0;
  } else if (depth >= height) {
    if (!is_open) {
      return 0.0;
    }
    return overflow_width;
  }
  return bottom_width + (cusp_width - bottom_width) * (depth / height);
}
double Cipolletti::wetted_perimeter(double depth) {
  if (depth <= 0.0) {
    return 0.0;
  } else if (depth >= height) {
    double a = (cusp_width - bottom_width) / 2.0;
    double side_length = sqrt(a * a + height * height);
    if (!is_open) {
      return bottom_width + cusp_width + 2.0 * side_length;
    }
    return bottom_width + 2.0 * side_length + 2.0 * (depth - height) +
           (overflow_width - cusp_width);
  }
  double t = (depth / height) * (cusp_width - bottom_width) + bottom_width;
  double a = (t - bottom_width) / 2.0;
  double side_length = sqrt(a * a + depth * depth);
  return 2.0 * side_length + bottom_width;
}
double Cipolletti::flow_area(double depth) {
  if (depth <= 0.0) {
    return 0.0;
  } else if (depth >= height) {
    if (!is_open) {
      return height * (cusp_width + bottom_width) / 2.0;
    }
    return (depth - height) * overflow_width +
           height * (cusp_width + bottom_width) / 2.0;
  }
  double t = (depth / height) * (cusp_width - bottom_width) + bottom_width;
  return depth * (t + bottom_width) / 2.0;
}

} // namespace hazen
