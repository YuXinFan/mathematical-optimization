#include <iostream>

#include "AllFunctors.hpp"
#include "DownhillSimplex.hpp"

int main( int argc, char** argv ) {

  PowellFunctor po;
  DownhillSimplex<PowellFunctor> optimizer(po);

  Eigen::Vector4d x_init;
  x_init << 1.0, 1.0, 1.0, 1.0;;
  Eigen::Vector4d x_opt;
  x_opt << 0.0, 0.0, 0.0, 0.0;
  Eigen::Vector4d x = x_init;
  int result = optimizer.minimize(x);

  double error = (x-x_opt).norm();
  if (error < 0.1)
    std::cout << "converged\n";
  else
    std::cout << "failed\n";

  return 0;
}