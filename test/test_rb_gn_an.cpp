#include <iostream>

#include "AllFunctors.hpp"
#include "GaussNewton.hpp"

int main( int argc, char** argv ) {

  RosenbrockFunctor rb;
  GaussNewton<RosenbrockFunctor> optimizer(rb);

  Eigen::Vector2d x_init;
  x_init << 1.0, 3.0;
  Eigen::Vector2d x_opt;
  x_opt << 1.0, 1.0;
  Eigen::Vector2d x = x_init;
  int result = optimizer.minimize(x);

  double error = (x-x_opt).norm();
  if (error < 0.1)
    std::cout << "converged\n";
  else
    std::cout << "failed\n";

  return 0;
}