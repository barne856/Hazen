#include "HydraulicShapes.hpp"
#include "HydraulicUtil.hpp"
#include <exception>
#include <iostream>
#define _USE_MATH_DEFINES
#include <math.h>

namespace hazen {

HydraulicShape::HydraulicShape(double height, bool is_open)
    : height(height), is_open(is_open) {}

double HydraulicShape::hydraulic_radius(double depth) {
  double P = wetted_perimeter(depth);
  if (P > 0.0) {
    return flow_area(depth) / wetted_perimeter(depth);
  }
  return 0.0;
}
double HydraulicShape::hydraulic_diameter(double depth) {
  return 4.0 * hydraulic_radius(depth);
}
double HydraulicShape::hydraulic_depth(double depth) {
  double T = top_width(depth);
  if (!is_open && depth >= height) {
    return std::numeric_limits<double>::infinity();
  }
  if (T > 0.0) {
    return flow_area(depth) / top_width(depth);
  }
  return 0.0f;
}
double HydraulicShape::get_max_depth() {
  if (is_open) {
    return std::numeric_limits<double>::infinity();
  }
  return height;
}
double HydraulicShape::froude(double Q, double depth) {
  if (depth <= 0.0) {
    return std::numeric_limits<double>::infinity();
  }
  double V = Q / flow_area(depth);
  double hm = hydraulic_depth(depth);
  return hazen::froude(V, hm);
}

Circle::Circle(double diameter) : HydraulicShape(diameter) {}
double Circle::top_width(double depth) {
  if (depth <= 0.0) {
    return 0.0;
  } else if (depth >= height) {
    return 0.0;
  }
  double T;
  try {
    T = 2.0 * sqrt((height * height) / 4.0 - pow((depth - height / 2.0), 2.0));
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    T = 0.0;
  }
  return T;
}

double Circle::wetted_perimeter(double depth) {
  if (depth >= height) {
    return M_PI * height;
  } else if (depth <= 0.0) {
    return 0.0;
  }
  double r = height / 2.0;
  double theta;
  try {
    theta = acos(1.0 - depth / r);
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    theta = M_PI;
  }
  return theta * height;
}

double Circle::flow_area(double depth) {
  if (depth >= height) {
    return M_PI * height * height / 4.0;
  } else if (depth <= 0.0) {
    return 0.0;
  }
  double r = height / 2.0;
  double theta;
  try {
    theta = acos(1.0 - depth / r);
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    theta = M_PI;
  }
  double A;
  try {
    A = r * r * theta - (r - depth) * sqrt(2.0 * r * depth - depth * depth);
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    A = M_PI * height * height / 4.0;
  }
  return A;
}

Rectangle::Rectangle(double width, double height, bool is_open)
    : length(width), HydraulicShape(height, is_open) {}

double Rectangle::top_width(double depth) {
  if (depth <= 0.0) {
    return 0.0;
  } else if (depth >= height && !is_open) {
    return 0.0;
  }
  return length;
}

double Rectangle::wetted_perimeter(double depth) {
  if (depth <= 0.0) {
    return 0.0;
  } else if (depth >= height && !is_open) {
    return 2.0 * length + 2.0 * height;
  }
  return length + 2.0 * depth;
}

double Rectangle::flow_area(double depth) {
  if (depth <= 0.0) {
    return 0.0;
  } else if (depth >= height && !is_open) {
    return length * height;
  }
  return length * depth;
}

VNotch::VNotch(double angle, double height, double overflow_length,
               bool is_open)
    : angle(angle), HydraulicShape(height, is_open),
      overflow_length(overflow_length) {}

double VNotch::top_width(double depth) {
  if (depth <= 0.0) {
    return 0.0;
  } else if (depth >= height) {
    if (!is_open) {
      return 0.0;
    }
    return overflow_length;
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
    return overflow_length - top_length + 2.0 * height / cos(rads / 2.0) +
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
    return overflow_length * (depth - height) + 0.5 * height * top_length;
  }
  double rads = angle * M_PI / 180.0;
  double top_length = 2.0 * depth * tan(rads / 2.0);
  return 0.5 * depth * top_length;
}

Cipolletti::Cipolletti(double top_length, double bottom_length, double height,
                       double overflow_length, bool is_open)
    : top_length(top_length), bottom_length(bottom_length),
      HydraulicShape(height, is_open), overflow_length(overflow_length) {}

double Cipolletti::top_width(double depth) {
  if (depth <= 0.0) {
    return 0.0;
  } else if (depth >= height) {
    if (!is_open) {
      return 0.0;
    }
    return overflow_length;
  }
  return bottom_length + (top_length - bottom_length) * (depth / height);
}

double Cipolletti::wetted_perimeter(double depth) {
  if (depth <= 0.0) {
    return 0.0;
  } else if (depth >= height) {
    double a = (top_length - bottom_length) / 2.0;
    double side_length = sqrt(a * a + height * height);
    if (!is_open) {
      return bottom_length + top_length + 2.0 * side_length;
    }
    return bottom_length + 2.0 * side_length + 2.0 * (depth - height) +
           (overflow_length - top_length);
  }
  double t = (depth / height) * (top_length - bottom_length) + bottom_length;
  double a = (t - bottom_length) / 2.0;
  double side_length = sqrt(a * a + depth * depth);
  return 2.0 * side_length + bottom_length;
}

double Cipolletti::flow_area(double depth) {
  if (depth <= 0.0) {
    return 0.0;
  } else if (depth >= height) {
    if (!is_open) {
      return height * (top_length + bottom_length) / 2.0;
    }
    return (depth - height) * overflow_length +
           height * (top_length + bottom_length) / 2.0;
  }
  double t = (depth / height) * (top_length - bottom_length) + bottom_length;
  return depth * (t + bottom_length) / 2.0;
}

} // namespace hazen
