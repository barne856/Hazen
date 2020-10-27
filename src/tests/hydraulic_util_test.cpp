#include "HeadLoss.hpp"
#include "HydraulicUtil.hpp"
#include <iostream>

double objective(std::vector<double> x) { return x[0] * x[1] - 1.0; }
std::vector<double> constraints(std::vector<double> x) { return {x[0] - x[1]}; }

using namespace hazen;

int main() {
  // std::vector<double> x0{0.0, 0.0};
  // auto result = hazen::constrained_solver_secant(x0, 0.000001, objective,
  // constraints); std::cout << "result = " << result[0] << ", " << result[1] <<
  // std::endl;

  HydraulicShape *channel_shape = new Rectangle(5.0, 10.0);
  double invert = 0.0;
  HydraulicShape *bar_screen_opening_shape = new Rectangle(0.5 / 12.0, 8.0);
  int n = 10;
  double Q = 10.0;
  double h = 4.0;

  double H =
      bar_screen_loss(channel_shape, invert, bar_screen_opening_shape, n, Q, h);

  std::cout << " The upstream Energy head of the bar screen is: " << H
            << std::endl;

  return 0;
}