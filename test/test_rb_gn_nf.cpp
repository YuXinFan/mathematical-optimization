#include <iostream>

#include "AllFunctors.hpp"
#include "NumDiff.hpp"
#include "GaussNewton.hpp"

int main( int argc, char** argv ) {

  RosenbrockFunctor rb;
  NumDiff<RosenbrockFunctor,Forward> num_rb(rb);
  GaussNewton< NumDiff<RosenbrockFunctor,Forward> > optimizer(num_rb);

  Eigen::Vector2d x_init;
  x_init << 0.0, 0.0;
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