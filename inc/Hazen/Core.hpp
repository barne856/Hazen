#ifndef HAZEN_CORE
#define HAZEN_CORE
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/LU>
#include <Eigen/QR>
#include <cmath>
#include <iostream>
#include <type_traits>
#include <utility>

namespace hazen {

// Implementation --------------------------------------------------------------
namespace Unit_impl {
using fraction = std::pair<int, int>;
constexpr int gcd(int a, int b) {
  if (a == 0)
    return b;
  return gcd(b % a, a);
}
constexpr fraction simplify(fraction frac) {
  int common_factor = gcd(frac.first, frac.second);
  fraction result = {frac.first / common_factor, frac.second / common_factor};
  if (result.second < 0) {
    result.first *= -1;
    result.second *= -1;
  }
  return result;
}
constexpr fraction add_fraction(fraction frac1, fraction frac2) {
  int den = gcd(frac1.second, frac2.second);
  den = frac1.second * frac2.second / den;
  int num =
      frac1.first * (den / frac1.second) + frac2.first * (den / frac2.second);
  return simplify({num, den});
}
constexpr fraction subtract_fraction(fraction frac1, fraction frac2) {
  frac2.first *= -1;
  return add_fraction(frac1, frac2);
}
} // namespace Unit_impl

// Unit Templates --------------------------------------------------------------
template <int L_num, int L_den, int T_num, int T_den, int M_num, int M_den,
          int Theta_num, int Theta_den>
struct Unit {
  static constexpr Unit_impl::fraction L = {L_num, L_den};
  static constexpr Unit_impl::fraction T = {T_num, T_den};
  static constexpr Unit_impl::fraction M = {M_num, M_den};
  static constexpr Unit_impl::fraction Theta = {Theta_num, Theta_den};
};

template <typename U1, typename U2> struct UPlus {
  using type = Unit<Unit_impl::add_fraction(U1::L, U2::L).first,
                    Unit_impl::add_fraction(U1::L, U2::L).second,
                    Unit_impl::add_fraction(U1::T, U2::T).first,
                    Unit_impl::add_fraction(U1::T, U2::T).second,
                    Unit_impl::add_fraction(U1::M, U2::M).first,
                    Unit_impl::add_fraction(U1::M, U2::M).second,
                    Unit_impl::add_fraction(U1::Theta, U2::Theta).first,
                    Unit_impl::add_fraction(U1::Theta, U2::Theta).second>;
};
template <typename U1, typename U2> struct UMinus {
  using type = Unit<Unit_impl::subtract_fraction(U1::L, U2::L).first,
                    Unit_impl::subtract_fraction(U1::L, U2::L).second,
                    Unit_impl::subtract_fraction(U1::T, U2::T).first,
                    Unit_impl::subtract_fraction(U1::T, U2::T).second,
                    Unit_impl::subtract_fraction(U1::M, U2::M).first,
                    Unit_impl::subtract_fraction(U1::M, U2::M).second,
                    Unit_impl::subtract_fraction(U1::Theta, U2::Theta).first,
                    Unit_impl::subtract_fraction(U1::Theta, U2::Theta).second>;
};
template <typename U, int divisor> struct UDivide {
  static constexpr Unit_impl::fraction L =
      Unit_impl::add_fraction({U::L.first, U::L.second *divisor}, {0, 1});
  static constexpr Unit_impl::fraction T =
      Unit_impl::add_fraction({U::T.first, U::T.second *divisor}, {0, 1});
  static constexpr Unit_impl::fraction M =
      Unit_impl::add_fraction({U::M.first, U::M.second *divisor}, {0, 1});
  static constexpr Unit_impl::fraction Theta = Unit_impl::add_fraction(
      {U::Theta.first, U::Theta.second *divisor}, {0, 1});
  using type = Unit<L.first, L.second, T.first, T.second, M.first, M.second,
                    Theta.first, Theta.second>;
};
template <typename U, int multiplier> struct UMultiply {
  static constexpr Unit_impl::fraction L =
      Unit_impl::add_fraction({U::L.first * multiplier, U::L.second}, {0, 1});
  static constexpr Unit_impl::fraction T =
      Unit_impl::add_fraction({U::T.first * multiplier, U::T.second}, {0, 1});
  static constexpr Unit_impl::fraction M =
      Unit_impl::add_fraction({U::M.first * multiplier, U::M.second}, {0, 1});
  static constexpr Unit_impl::fraction Theta = Unit_impl::add_fraction(
      {U::Theta.first * multiplier, U::Theta.second}, {0, 1});
  using type = Unit<L.first, L.second, T.first, T.second, M.first, M.second,
                    Theta.first, Theta.second>;
};

template <typename U1, typename U2>
using Unit_Plus = typename UPlus<U1, U2>::type;
template <typename U1, typename U2>
using Unit_Minus = typename UMinus<U1, U2>::type;
template <typename U, int divisor>
using Unit_Divide = typename UDivide<U, divisor>::type;
template <typename U, int multiplier>
using Unit_Multiply = typename UMultiply<U, multiplier>::type;

// Unit Definitions ------------------------------------------------------------
using Dimensionless_Unit = Unit<0, 1, 0, 1, 0, 1, 0, 1>;
using Length_Unit = Unit<1, 1, 0, 1, 0, 1, 0, 1>;
using Time_Unit = Unit<0, 1, 1, 1, 0, 1, 0, 1>;
using Mass_Unit = Unit<0, 1, 0, 1, 1, 1, 0, 1>;
using Temperature_Unit = Unit<0, 1, 0, 1, 0, 1, 1, 1>;
using Velocity_Unit = Unit_Minus<Length_Unit, Time_Unit>;
using Acceleration_Unit = Unit_Minus<Velocity_Unit, Time_Unit>;
using Area_Unit = Unit_Plus<Length_Unit, Length_Unit>;
using Volume_Unit = Unit_Plus<Area_Unit, Length_Unit>;
using Density_Unit = Unit_Minus<Mass_Unit, Volume_Unit>;
using Force_Unit = Unit_Plus<Mass_Unit, Acceleration_Unit>;
using Pressure_Unit = Unit_Minus<Force_Unit, Area_Unit>;
using Dynamic_Viscosity_Unit = Unit_Plus<Pressure_Unit, Time_Unit>;
using Flow_Unit = Unit_Minus<Volume_Unit, Time_Unit>;

// Scalar Templates ------------------------------------------------------------
template <typename U> class Scalar {
public:
  constexpr Scalar() : val{} {}
  explicit constexpr Scalar(double d) : val{d} {}
  Scalar<U> &operator+=(const Scalar<U> &x) {
    val += x.val;
    return *this;
  }
  Scalar<U> &operator-=(const Scalar<U> &x) {
    val -= x.val;
    return *this;
  }
  friend std::ostream &operator<<(std::ostream &os, const Scalar<U> &s) {
    os << s.val;
    return os;
  }

  double val;
};

// Binary Arithmetic Operations
template <typename U> inline Scalar<U> operator+(Scalar<U> x, Scalar<U> y) {
  return Scalar<U>{x.val + y.val};
}
template <typename U> inline Scalar<U> operator-(Scalar<U> x, Scalar<U> y) {
  return Scalar<U>{x.val - y.val};
}
template <typename U1, typename U2>
inline Scalar<Unit_Plus<U1, U2>> operator*(Scalar<U1> x, Scalar<U2> y) {
  return Scalar<Unit_Plus<U1, U2>>{x.val * y.val};
}
template <typename U1, typename U2>
inline Scalar<Unit_Minus<U1, U2>> operator/(Scalar<U1> x, Scalar<U2> y) {
  return Scalar<Unit_Minus<U1, U2>>{x.val / y.val};
}

// Use doubles as dimensionless scalars
template <typename U> inline Scalar<U> operator*(Scalar<U> x, double y) {
  return Scalar<U>{x.val * y};
}
template <typename U> inline Scalar<U> operator*(double x, Scalar<U> y) {
  return Scalar<U>{x * y.val};
}
template <typename U> inline Scalar<U> operator/(Scalar<U> x, double y) {
  return Scalar<U>{x.val / y};
}
template <typename U>
inline Scalar<Unit_Minus<Unit<0, 1, 0, 1, 0, 1, 0, 1>, U>>
operator/(double x, Scalar<U> y) {
  return Scalar<Unit_Minus<Unit<0, 1, 0, 1, 0, 1, 0, 1>, U>>{x / y.val};
}

// Binary Comparison Operations
template <typename U> inline bool operator==(Scalar<U> x, Scalar<U> y) {
  return x.val == y.val;
}
template <typename U> inline bool operator!=(Scalar<U> x, Scalar<U> y) {
  return x.val != y.val;
}
template <typename U> inline bool operator>=(Scalar<U> x, Scalar<U> y) {
  return x.val >= y.val;
}
template <typename U> inline bool operator<=(Scalar<U> x, Scalar<U> y) {
  return x.val <= y.val;
}
template <typename U> inline bool operator>(Scalar<U> x, Scalar<U> y) {
  return x.val > y.val;
}
template <typename U> inline bool operator<(Scalar<U> x, Scalar<U> y) {
  return x.val < y.val;
}

// Scalar Functions ------------------------------------------------------------
template <int pow_num, int pow_den, typename U>
constexpr inline Scalar<Unit_Divide<Unit_Multiply<U, pow_num>, pow_den>>
power(Scalar<U> x) {
  return Scalar<Unit_Divide<Unit_Multiply<U, pow_num>, pow_den>>(std::pow(
      x.val, static_cast<double>(pow_num) / static_cast<double>(pow_den)));
}

template <typename U>
constexpr inline Scalar<Unit_Divide<U, 2>> sqrt(Scalar<U> x) {
  return Scalar<Unit_Divide<U, 2>>(std::sqrt(x.val));
}
template <typename U> constexpr inline Scalar<U> abs(Scalar<U> x) {
  return Scalar<U>(std::abs(x.val));
}
template <typename U>
constexpr inline Scalar<Unit_Multiply<U, 2>> abs2(Scalar<U> x) {
  return x * x;
}

// Scalar Result Types ---------------------------------------------------------
// Product type of two scalars
template <typename T1, typename T2>
using product_type = decltype(std::declval<T1 &>() * std::declval<T2 &>());
// Quotient type of two scalars
template <typename T1, typename T2>
using quotient_type = decltype(std::declval<T1 &>() / std::declval<T2 &>());
// Root of a scalar
template <typename T> using root_type = decltype(sqrt(std::declval<T>()));

// Scalar Definitions ----------------------------------------------------------
using Dimensionless = Scalar<Dimensionless_Unit>;
using Angle = Scalar<Dimensionless_Unit>;
using Length = Scalar<Length_Unit>;
using Time = Scalar<Time_Unit>;
using Mass = Scalar<Mass_Unit>;
using Temperature = Scalar<Temperature_Unit>;
using Velocity = Scalar<Velocity_Unit>;
using Acceleration = Scalar<Acceleration_Unit>;
using Area = Scalar<Area_Unit>;
using Volume = Scalar<Volume_Unit>;
using Density = Scalar<Density_Unit>;
using Force = Scalar<Force_Unit>;
using Pressure = Scalar<Pressure_Unit>;
using Dynamic_Viscosity = Scalar<Dynamic_Viscosity_Unit>;
using Flow = Scalar<Flow_Unit>;

template <> class Scalar<Dimensionless_Unit> {
public:
  constexpr Scalar() : val{} {}
  explicit constexpr Scalar(double d) : val{d} {}
  static constexpr Angle Radians(double d) { return Angle(d); }
  static constexpr Angle Degrees(double d) { return Angle(d * M_PI / 180.0); }
  static Angle Slope(double rise, double run) {
    return Angle(std::atan2(rise, run));
  }
  static Angle Slope(double slope) { return Angle(std::atan(slope)); }
  static Angle Percent(double d) { return Angle(std::atan(d * 0.01)); }
  constexpr double as_radians() const { return val; }
  constexpr double as_degrees() const { return val * 180.0 / M_PI; }
  double as_slope() const { return std::tan(val); }
  double as_percent() const { return std::tan(val) * 100.0; }

  Scalar &operator+=(const Scalar &x) {
    val += x.val;
    return *this;
  }
  Scalar &operator-=(const Scalar &x) {
    val -= x.val;
    return *this;
  }
  friend std::ostream &operator<<(std::ostream &os, const Scalar &s) {
    os << s.val;
    return os;
  }

  double val;
};
template <> class Scalar<Length_Unit> {
public:
  constexpr Scalar() : val{} {}
  explicit constexpr Scalar(double d) : val{d} {}
  static constexpr Length Feet(double d) { return Length(d * 0.3048); }
  static constexpr Length Yards(double d) { return Length(d * 0.3048 * 3); }
  static constexpr Length Inches(double d) { return Length(d * 0.3048 / 12); }
  static constexpr Length Meters(double d) { return Length(d); }
  constexpr double as_feet() const { return val / 0.3048; }
  constexpr double as_yards() const { return val / 0.3048 / 3; }
  constexpr double as_inches() const { return val / 0.3048 * 12; }
  constexpr double as_meters() const { return val; }

  Scalar &operator+=(const Scalar &x) {
    val += x.val;
    return *this;
  }
  Scalar &operator-=(const Scalar &x) {
    val -= x.val;
    return *this;
  }
  friend std::ostream &operator<<(std::ostream &os, const Scalar &s) {
    os << s.val;
    return os;
  }

  double val;
};
template <> class Scalar<Time_Unit> {
public:
  constexpr Scalar() : val{} {}
  explicit constexpr Scalar(double d) : val{d} {}
  static constexpr Time Seconds(double d) { return Time(d); }
  constexpr double as_seconds() const { return val; }
  constexpr double as_minutes() const { return val / 60.0; }
  constexpr double as_hours() const { return val / 60.0 / 60.0; }
  constexpr double as_days() const { return val / 60.0 / 60.0 / 24.0; }

  Scalar &operator+=(const Scalar &x) {
    val += x.val;
    return *this;
  }
  Scalar &operator-=(const Scalar &x) {
    val -= x.val;
    return *this;
  }
  friend std::ostream &operator<<(std::ostream &os, const Scalar &s) {
    os << s.val;
    return os;
  }

  double val;
};
template <> class Scalar<Mass_Unit> {
public:
  constexpr Scalar() : val{} {}
  explicit constexpr Scalar(double d) : val{d} {}
  static constexpr Mass Kilograms(double d) { return Mass(d); }
  static constexpr Mass Slug(double d) { return Mass(d * 14.593903); }
  constexpr double as_kilogram() const { return val; }
  constexpr double as_slug() const { return val / 14.593903; }

  Scalar &operator+=(const Scalar &x) {
    val += x.val;
    return *this;
  }
  Scalar &operator-=(const Scalar &x) {
    val -= x.val;
    return *this;
  }
  friend std::ostream &operator<<(std::ostream &os, const Scalar &s) {
    os << s.val;
    return os;
  }

  double val;
};
template <> class Scalar<Temperature_Unit> {
public:
  constexpr Scalar() : val{} {}
  explicit constexpr Scalar(double d) : val{d} {}
  static constexpr Temperature Kelvin(double d) { return Temperature(d); }
  static constexpr Temperature Celcius(double d) {
    return Temperature(d + 273.15);
  }
  static constexpr Temperature Fahrenheit(double d) {
    return Temperature((d - 32) * 5 / 9 + 237.15);
  }
  constexpr double as_kelvin(double d) const { return val; }
  constexpr double as_celcius(double d) const { return val - 273.15; }
  constexpr double as_fahrenheit(double d) const {
    return (val - 273.15) * 9 / 5 + 32;
  }

  Scalar &operator+=(const Scalar &x) {
    val += x.val;
    return *this;
  }
  Scalar &operator-=(const Scalar &x) {
    val -= x.val;
    return *this;
  }
  friend std::ostream &operator<<(std::ostream &os, const Scalar &s) {
    os << s.val;
    return os;
  }

  double val;
};
template <> class Scalar<Velocity_Unit> {
public:
  constexpr Scalar() : val{} {}
  explicit constexpr Scalar(double d) : val{} {}
  static constexpr Velocity MPS(double d) { return Velocity(d); }
  static constexpr Velocity FPS(double d) { return Velocity(d * 0.3048); }
  constexpr double as_meters_per_second() const { return val; }
  constexpr double as_feet_per_second() const { return val / 0.3048; }

  Scalar &operator+=(const Scalar &x) {
    val += x.val;
    return *this;
  }
  Scalar &operator-=(const Scalar &x) {
    val -= x.val;
    return *this;
  }
  friend std::ostream &operator<<(std::ostream &os, const Scalar &s) {
    os << s.val;
    return os;
  }

  double val;
};
template <> class Scalar<Acceleration_Unit> {
public:
  constexpr Scalar() : val{} {}
  explicit constexpr Scalar(double d) : val{d} {}
  static constexpr Acceleration MPS2(double d) { return Acceleration(d); }
  static constexpr Acceleration FPS2(double d) {
    return Acceleration(d * 0.3048);
  }
  constexpr double as_meters_per_second_squared() const { return val; }
  constexpr double as_feet_per_second_squared() const { return val / 0.3048; }

  Scalar &operator+=(const Scalar &x) {
    val += x.val;
    return *this;
  }
  Scalar &operator-=(const Scalar &x) {
    val -= x.val;
    return *this;
  }
  friend std::ostream &operator<<(std::ostream &os, const Scalar &s) {
    os << s.val;
    return os;
  }

  double val;
};
template <> class Scalar<Area_Unit> {
public:
  constexpr Scalar() : val{} {}
  explicit constexpr Scalar(double d) : val{d} {}
  static constexpr Area Square_Meters(double d) { return Area(d); }
  static constexpr Area Square_Feet(double d) {
    return Area(d * 0.3048 * 0.3048);
  }
  static constexpr Area Square_Inches(double d) {
    return Area(d * 0.3048 * 0.3048 / 12 / 12);
  }
  static constexpr Area Square_Yards(double d) {
    return Area(d * 0.3048 * 0.3048 * 3 * 3);
  }
  static constexpr Area Acres(double d) { return Area(d * 4046.8564224); }
  constexpr double as_square_meters() const { return val; }
  constexpr double as_square_feet() const { return val / 0.3048 / 0.3048; }
  constexpr double as_square_inches() const {
    return 12 * 12 * val / 0.3048 / 0.3048;
  }
  constexpr double as_square_yards() const {
    return val / 0.3048 / 0.3048 / 3 / 3;
  }
  constexpr double as_acres() const { return val / 4046.8564224; }

  Scalar &operator+=(const Scalar &x) {
    val += x.val;
    return *this;
  }
  Scalar &operator-=(const Scalar &x) {
    val -= x.val;
    return *this;
  }
  friend std::ostream &operator<<(std::ostream &os, const Scalar &s) {
    os << s.val;
    return os;
  }

  double val;
};
template <> class Scalar<Volume_Unit> {
public:
  constexpr Scalar() : val{} {}
  explicit constexpr Scalar(double d) : val{d} {}
  static constexpr Volume Cubic_Meters(double d) { return Volume(d); }
  static constexpr Volume Cubic_Feet(double d) {
    return Volume(d * 0.3048 * 0.3048 * 0.3048);
  }
  static constexpr Volume Gallons(double d) { return Volume(d * 0.0037854118); }
  static constexpr Volume Liters(double d) { return Volume(d * 0.001); }
  constexpr double as_cubic_meters() const { return val; }
  constexpr double as_cubic_feet() const {
    return val / 0.3048 / 0.3048 / 0.3048;
  }
  constexpr double as_gallons() const { return val / 0.0037854118; }
  constexpr double as_liters() const { return val / 0.001; }

  Scalar &operator+=(const Scalar &x) {
    val += x.val;
    return *this;
  }
  Scalar &operator-=(const Scalar &x) {
    val -= x.val;
    return *this;
  }
  friend std::ostream &operator<<(std::ostream &os, const Scalar &s) {
    os << s.val;
    return os;
  }

  double val;
};
template <> class Scalar<Density_Unit> {
public:
  constexpr Scalar() : val{} {}
  explicit constexpr Scalar(double d) : val{d} {}
  static constexpr Density Kilograms_per_Cubic_Meter(double d) {
    return Density(d);
  }
  static constexpr Density Slugs_per_Cubic_Foot(double d) {
    return Density(d * 515.3788184);
  }
  constexpr double as_kilograms_per_cubic_meter() { return val; }
  constexpr double as_slugs_per_cubic_foot() { return val / 515.3788184; }

  Scalar &operator+=(const Scalar &x) {
    val += x.val;
    return *this;
  }
  Scalar &operator-=(const Scalar &x) {
    val -= x.val;
    return *this;
  }
  friend std::ostream &operator<<(std::ostream &os, const Scalar &s) {
    os << s.val;
    return os;
  }

  double val;
};
template <> class Scalar<Force_Unit> {
public:
  constexpr Scalar() : val{} {}
  explicit constexpr Scalar(double d) : val{d} {}
  static constexpr Force Newtons(double d) { return Force(d); }
  static constexpr Force Pounds(double d) { return Force(d * 4.4482216); }
  constexpr double as_newtons(double d) const { return val; }
  constexpr double as_pounds(double d) const { return val / 4.4482216; }

  Scalar &operator+=(const Scalar &x) {
    val += x.val;
    return *this;
  }
  Scalar &operator-=(const Scalar &x) {
    val -= x.val;
    return *this;
  }
  friend std::ostream &operator<<(std::ostream &os, const Scalar &s) {
    os << s.val;
    return os;
  }

  double val;
};
template <> class Scalar<Pressure_Unit> {
public:
  constexpr Scalar() : val{} {}
  explicit constexpr Scalar(double d) : val{d} {}
  static constexpr Pressure Pascals(double d) { return Pressure(d); }
  static constexpr Pressure PSF(double d) {
    return Pressure(d * 47.880258888889);
  }
  static constexpr Pressure PSI(double d) { return Pressure(d * 6894.7572932); }
  constexpr double as_pascals() const { return val; }
  constexpr double as_psf() const { return val / 47.880258888889; }
  constexpr double as_psi() const { return val / 6894.7572932; }

  Scalar &operator+=(const Scalar &x) {
    val += x.val;
    return *this;
  }
  Scalar &operator-=(const Scalar &x) {
    val -= x.val;
    return *this;
  }
  friend std::ostream &operator<<(std::ostream &os, const Scalar &s) {
    os << s.val;
    return os;
  }

  double val;
};
template <> class Scalar<Dynamic_Viscosity_Unit> {
public:
  constexpr Scalar() : val{} {}
  explicit constexpr Scalar(double d) : val{d} {}
  static constexpr Dynamic_Viscosity Pascal_Seconds(double d) {
    return Dynamic_Viscosity(d);
  }
  static constexpr Dynamic_Viscosity PSF_Seconds(double d) {
    return Dynamic_Viscosity(d * 47.880258888889);
  }
  constexpr double as_pascal_seconds() const { return val; }
  constexpr double as_psf_seconds() const { return val / 47.880258888889; }

  Scalar &operator+=(const Scalar &x) {
    val += x.val;
    return *this;
  }
  Scalar &operator-=(const Scalar &x) {
    val -= x.val;
    return *this;
  }
  friend std::ostream &operator<<(std::ostream &os, const Scalar &s) {
    os << s.val;
    return os;
  }

  double val;
};
template <> class Scalar<Flow_Unit> {
public:
  constexpr Scalar() : val{} {}
  explicit constexpr Scalar(double d) : val{d} {}
  static constexpr Flow CMS(double d) { return Flow(d); }
  static constexpr Flow LPS(double d) { return Flow(d * 0.001); }
  static constexpr Flow CFS(double d) {
    return Flow(d * 0.3048 * 0.3048 * 0.3048);
  }
  static constexpr Flow GPM(double d) { return Flow(d * 0.0037854118 / 60); }
  static constexpr Flow GPD(double d) {
    return Flow(d * 0.0037854118 / 60 / 60 / 24);
  }
  static constexpr Flow MGD(double d) {
    return Flow(d * 0.0037854118 / 60 / 60 / 24 * 1000000);
  }
  constexpr double as_cubic_meters_per_second() const { return val; }
  constexpr double as_liters_per_second() const { return val / 0.001; }
  constexpr double as_cubic_feet_per_second() const {
    return val / 0.3048 / 0.3048 / 0.3048;
  }
  constexpr double as_gallons_per_minute() const {
    return val * 60 / 0.0037854118;
  }
  constexpr double as_gallons_per_day() const {
    return val * 60 * 60 * 24 / 0.0037854118;
  }
  constexpr double as_million_gallons_per_day() const {
    return val * 60 * 60 * 24 / 0.0037854118 / 1000000;
  }

  Scalar &operator+=(const Scalar &x) {
    val += x.val;
    return *this;
  }
  Scalar &operator-=(const Scalar &x) {
    val -= x.val;
    return *this;
  }
  friend std::ostream &operator<<(std::ostream &os, const Scalar &s) {
    os << s.val;
    return os;
  }

  double val;
};

// Unit Literals ---------------------------------------------------------------
// Meters
constexpr Length operator""_m(long double d) {
  return Length::Meters(static_cast<double>(d));
}
// Yard
constexpr Length operator""_yd(long double d) {
  return Length::Yards(static_cast<double>(d));
}
// Feet
constexpr Length operator""_ft(long double d) {
  return Length::Feet(static_cast<double>(d));
}
// Inches
constexpr Length operator""_in(long double d) {
  return Length::Inches(static_cast<double>(d));
}
// Kilograms
constexpr Mass operator""_kg(long double d) {
  return Mass::Kilograms(static_cast<double>(d));
}
// Slug
constexpr Mass operator""_slug(long double d) {
  return Mass::Slug(static_cast<double>(d));
}
// Seconds
constexpr Time operator""_s(long double d) {
  return Time::Seconds(static_cast<double>(d));
}
// Kelvin
constexpr Temperature operator""_K(long double d) {
  return Temperature::Kelvin(static_cast<double>(d));
}
// Celcius
constexpr Temperature operator""_C(long double d) {
  return Temperature::Celcius(static_cast<double>(d));
}
// Fahrenheit
constexpr Temperature operator""_F(long double d) {
  return Temperature::Fahrenheit(static_cast<double>(d));
}
// Cubic Meters
constexpr Volume operator""_m3(long double d) {
  return Volume::Cubic_Meters(static_cast<double>(d));
}
// Cubic feet
constexpr Volume operator""_ft3(long double d) {
  return Volume::Cubic_Feet(static_cast<double>(d));
}
// Liters
constexpr Volume operator""_l(long double d) {
  return Volume::Liters(static_cast<double>(d));
}
// Gallons
constexpr Volume operator""_gal(long double d) {
  return Volume::Gallons(static_cast<double>(d));
}
// Cubic Meters per Second
constexpr Flow operator""_cms(long double d) {
  return Flow::CMS(static_cast<double>(d));
}
// Cubic Feet per Second
constexpr Flow operator""_cfs(long double d) {
  return Flow::CFS(static_cast<double>(d));
}
// Gallons per Minute
constexpr Flow operator""_gpm(long double d) {
  return Flow::GPM(static_cast<double>(d));
}
// Gallons per Day
constexpr Flow operator""_gpd(long double d) {
  return Flow::GPD(static_cast<double>(d));
}
// Million Gallons per Day
constexpr Flow operator""_mgd(long double d) {
  return Flow::MGD(static_cast<double>(d));
}
// Meters per Second
constexpr Velocity operator""_mps(long double d) {
  return Velocity::MPS(static_cast<double>(d));
}
// Feet per Second
constexpr Velocity operator""_fps(long double d) {
  return Velocity::FPS(static_cast<double>(d));
}
// Meters per Second Squared
constexpr Acceleration operator""_mps2(long double d) {
  return Acceleration::MPS2(static_cast<double>(d));
}
// Feet per Second Squared
constexpr Acceleration operator""_fps2(long double d) {
  return Acceleration::FPS2(static_cast<double>(d));
}
// Kilograms per Cubic Meter
constexpr Density operator""_kgm3(long double d) {
  return Density::Kilograms_per_Cubic_Meter(static_cast<double>(d));
}
// Slug per Cubic Feet
constexpr Density operator""_slugcf(long double d) {
  return Density::Slugs_per_Cubic_Foot(static_cast<double>(d));
}
// Pascal Seconds
constexpr Dynamic_Viscosity operator""_pas(long double d) {
  return Dynamic_Viscosity::Pascal_Seconds(static_cast<double>(d));
}
// Pounds per Square Foot Second
constexpr Dynamic_Viscosity operator""_psfs(long double d) {
  return Dynamic_Viscosity::PSF_Seconds(static_cast<double>(d));
}
// Newtons
constexpr Force operator""_n(long double d) {
  return Force::Newtons(static_cast<double>(d));
}
// Pounds (Force)
constexpr Force operator""_lbf(long double d) {
  return Force::Pounds(static_cast<double>(d));
}
// Pounds (Force) per Cubic Foot
constexpr Scalar<Unit_Minus<Force_Unit, Volume_Unit>>
operator""_pcf(long double d) {
  return Scalar<Unit_Minus<Force_Unit, Volume_Unit>>{
      static_cast<double>(d * 157.0874606377)};
}
// Radians
constexpr Angle operator""_rads(long double d) {
  return Angle::Radians(static_cast<double>(d));
}
// Degrees
constexpr Angle operator""_degrees(long double d) {
  return Angle::Degrees(static_cast<double>(d));
}
// Dimensionless
constexpr Dimensionless operator""_pure(long double d) {
  return Dimensionless(static_cast<double>(d));
}

template <typename T, int rows, int cols> class Matrix {
public:
  using value_type = T;

  Matrix() = default;
  Matrix(const Matrix &) = default;
  Matrix(Matrix &&) = default;
  Matrix &operator=(const Matrix &) = default;
  Matrix &operator=(Matrix &&) = default;

  Matrix(size_t m, size_t n)
      : elems{Eigen::Matrix<double, rows, cols>::Zero(m, n)} {}
  Matrix(size_t n) : elems{Eigen::Matrix<double, rows, cols>::Zero(n)} {}
  Matrix(Eigen::Matrix<double, rows, cols> m) : elems{m} {}
  Matrix &operator=(const T &s) {
    elems.setConstant(s.val);
    return *this;
  }
  Matrix &operator=(const Eigen::Matrix<double, rows, cols> &m) {
    elems = m;
    return *this;
  }
  Matrix &operator+=(const Matrix &m) {
    elems += m.elems;
    return *this;
  }
  Matrix &operator-=(const Matrix &m) {
    elems -= m.elems;
    return *this;
  }
  Matrix &operator-() {
    elems = -elems;
    return *this;
  }
  friend std::ostream &operator<<(std::ostream &os, const Matrix &m) {
    os << m.elems;
    return os;
  }

  operator T &() { return elems(0, 0); }
  operator T() const { return elems(0, 0); }
  size_t size() const { return elems.size(); }
  size_t n_rows() const { return elems.rows(); }
  size_t n_cols() const { return elems.cols(); }

  T operator()(size_t i, size_t j = 0) const { return T(elems(i, j)); }

  static Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Constant(int m, int n,
                                                            const T &s) {
    return Matrix<T, Eigen::Dynamic, Eigen::Dynamic>{
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>::Constant(m, n,
                                                                        s.val)};
  }
  static Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Random(int m, int n) {
    return Matrix<T, Eigen::Dynamic, Eigen::Dynamic>{
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>::Random(m, n)};
  }
  static Matrix<T, Eigen::Dynamic, 1> Random(int n) {
    return Matrix<T, Eigen::Dynamic, 1>{
        Eigen::Matrix<double, Eigen::Dynamic, 1>::Random(n)};
  }

  Eigen::Matrix<double, rows, cols> elems;
};

// Binary Arithmetic Operations
template <typename T, int rows, int cols>
inline Matrix<T, rows, cols> operator+(const Matrix<T, rows, cols> &x,
                                       const Matrix<T, rows, cols> &y) {
  return Matrix<T, rows, cols>{x.elems + y.elems};
}
template <typename T, int rows, int cols>
inline Matrix<T, rows, cols> operator-(const Matrix<T, rows, cols> &x,
                                       const Matrix<T, rows, cols> &y) {
  return Matrix<T, rows, cols>{x.elems - y.elems};
}
template <typename T1, typename T2, int m, int n, int o>
inline Matrix<product_type<T1, T2>, m, o> operator*(const Matrix<T1, m, n> &x,
                                                    const Matrix<T2, n, o> &y) {
  return Matrix<product_type<T1, T2>, m, o>{x.elems * y.elems};
}

// Use doubles as dimensionless scalars
template <typename T, int rows, int cols>
inline Matrix<T, rows, cols> operator*(const Matrix<T, rows, cols> &x,
                                       double y) {
  return Matrix<T, rows, cols>{x.elems * y};
}
template <typename T, int rows, int cols>
inline Matrix<T, rows, cols> operator*(double x,
                                       const Matrix<T, rows, cols> &y) {
  return Matrix<T, rows, cols>{x * y.elems};
}
template <typename T, int rows, int cols>
inline Matrix<T, rows, cols> operator/(const Matrix<T, rows, cols> &x,
                                       double y) {
  return Matrix<T, rows, cols>{x.elems / y};
}
template <typename T, int rows, int cols>
inline Matrix<quotient_type<Dimensionless, T>, rows, cols>
operator/(double x, const Matrix<T, rows, cols> &y) {
  return Matrix<quotient_type<Dimensionless, T>, rows, cols>{x / y.elems};
}

// Scaling
template <typename T, typename U, int rows, int cols>
inline Matrix<product_type<T, Scalar<U>>, rows, cols>
operator*(const Matrix<T, rows, cols> &x, Scalar<U> y) {
  return Matrix<product_type<T, Scalar<U>>, rows, cols>{x.elems * y.val};
}
template <typename T, typename U, int rows, int cols>
inline Matrix<product_type<T, Scalar<U>>, rows, cols>
operator*(Scalar<U> x, const Matrix<T, rows, cols> &y) {
  return Matrix<product_type<T, Scalar<U>>, rows, cols>{x.val * y.elems};
}
template <typename T, typename U, int rows, int cols>
inline Matrix<quotient_type<T, Scalar<U>>, rows, cols>
operator/(const Matrix<T, rows, cols> &x, Scalar<U> y) {
  return Matrix<quotient_type<T, Scalar<U>>, rows, cols>{x.elems / y.val};
}
template <typename T, typename U, int rows, int cols>
inline Matrix<quotient_type<Scalar<U>, T>, rows, cols>
operator/(Scalar<U> x, const Matrix<T, rows, cols> &y) {
  return Matrix<quotient_type<Scalar<U>, T>, rows, cols>{x.val / y.elems};
}

// Binary Comparison Operations
template <typename T, int rows, int cols>
inline bool operator==(const Matrix<T, rows, cols> &x,
                       const Matrix<T, rows, cols> &y) {
  return x.elems == y.elems;
}
template <typename T, int rows, int cols>
inline bool operator!=(const Matrix<T, rows, cols> &x,
                       const Matrix<T, rows, cols> &y) {
  return x.elems != y.elems;
}
template <typename T, int rows, int cols>
inline bool operator>=(const Matrix<T, rows, cols> &x,
                       const Matrix<T, rows, cols> &y) {
  return x.elems >= y.elems;
}
template <typename T, int rows, int cols>
inline bool operator<=(const Matrix<T, rows, cols> &x,
                       const Matrix<T, rows, cols> &y) {
  return x.elems <= y.elems;
}
template <typename T, int rows, int cols>
inline bool operator>(const Matrix<T, rows, cols> &x,
                      const Matrix<T, rows, cols> &y) {
  return x.elems > y.elems;
}
template <typename T, int rows, int cols>
inline bool operator<(const Matrix<T, rows, cols> &x,
                      const Matrix<T, rows, cols> &y) {
  return x.elems < y.elems;
}

// Matrix Functions ------------------------------------------------------------
template <typename T, int rows, int cols>
T norm(const Matrix<T, rows, cols> &m) {
  return T{m.elems.norm()};
}
template <typename T, int rows, int cols>
Matrix<T, rows, cols> abs(const Matrix<T, rows, cols> &m) {
  return Matrix<T, rows, cols>{m.elems.cwiseAbs()};
}

template <typename T, int rows, int cols>
product_type<T, T> squared_norm(const Matrix<T, rows, cols> &m) {
  return product_type<T, T>{m.elems.squaredNorm()};
}
template <typename T, int rows, int cols>
Matrix<T, cols, rows> transpose(const Matrix<T, rows, cols> &m) {
  return m.elems.transpose();
}
template <typename T1, typename T2, int rows>
product_type<T1, T2> dot_product(const Matrix<T1, rows, 1> &u,
                                 const Matrix<T2, rows, 1> &v) {
  return product_type<T1, T2>(u.elems.dot(v.elems));
}
template <typename T1, typename T2, int n, int m>
Matrix<quotient_type<T2, T1>, m, 1> solve(const Matrix<T1, n, m> &A,
                                          const Matrix<T2, n, 1> &b) {
  return Matrix<quotient_type<T2, T1>, m, 1>{A.elems.lu().solve(b.elems)};
}
template <typename T1, typename T2, int n, int m>
Matrix<quotient_type<T2, T1>, m, 1>
solve_least_squares(const Matrix<T1, n, m> &A, const Matrix<T2, n, 1> &b) {
  // return Matrix<quotient_type<T2, T1>, m, 1>{
  //     A.elems.householderQr().solve(b.elems)};
  return Matrix<quotient_type<T2, T1>, m, 1>{
      A.elems.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b.elems)};
}

// Matrix and Vector aliases
template <typename T> using Mat = Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
template <typename T> using Vec = Matrix<T, Eigen::Dynamic, 1>;

} // namespace hazen

#endif