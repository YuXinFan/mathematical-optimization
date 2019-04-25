#include <iostream>

#include "AllFunctors.hpp"
#include "NumDiff.hpp"
#include "GaussNewton.hpp"

int main( int argc, char** argv ) {

  std::list<Eigen::VectorXd> ps;
  for( int i = 0; i < 1000; i++ )
    ps.push_back( Eigen::VectorXd::Random(3) );

  WeberFunctor we(ps);
  NumDiff<WeberFunctor,Forward> num_we(we);
  GaussNewton< NumDiff<WeberFunctor,Forward> > optimizer(num_we);

  Eigen::VectorXd x_init(3);
  x_init << 1.0, 1.0, 1.0;
  Eigen::VectorXd x_opt(3);
  x_opt << 0.0, 0.0, 0.0;
  Eigen::VectorXd x = x_init;
  int result = optimizer.minimize(x);

  double error = (x-x_opt).norm();
  if (error < 0.1)
    std::cout << "converged\n";
  else
    std::cout << "failed\n";

  return 0;
}