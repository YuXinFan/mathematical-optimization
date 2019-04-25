#include <iostream>

#include "AllFunctors.hpp"
#include "NumDiff.hpp"
#include "GaussNewton.hpp"

int main( int argc, char** argv ) {

  PowellFunctor po;
  NumDiff<PowellFunctor,Central> num_po(po);
  GaussNewton< NumDiff<PowellFunctor,Central> > optimizer(num_po);

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