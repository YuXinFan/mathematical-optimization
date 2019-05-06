#include "functions.hpp"

#include <cmath>

Eigen::Vector2d rosenbrock(const Eigen::Vector2d &x) {
  double f1 = 10*(x(1)-x(0)*x(0));
  double f2 = 1-x(0);
  return Eigen::Vector2d(f1, f2);
}

Eigen::Matrix2d rosenbrockJ(const Eigen::Vector2d &x) {
  Eigen::Matrix2d jac;
  jac << -20*x(0), 10,
         -1, 0;
  return jac; 
}

Eigen::Vector4d powell(const Eigen::Vector4d &x) {
  double f1 = x(0) + 10 * x(1);
  double f2 = sqrt(5)*(x(2)-x(3));
  double f3 = (x(1)-2*x(2))*(x(1)-2*x(2));
  double f4 = sqrt(10)*((x(0)-x(3))*(x(0)-x(3)));
  return Eigen::Vector4d(f1, f2, f3, f4);
}

Eigen::Matrix4d powellJ(const Eigen::Vector4d &x) {
  Eigen::Matrix4d jac;
  jac << 1, 10, 0, 0, 
         0, 0, sqrt(5), -sqrt(5),
         0, 2*x(1)-4*x(2), -4*x(1)+8*x(2), 0,
         sqrt(10)*2*(x(0)-x(3)), 0, 0, sqrt(10)*2*(x(3)-x(0));
  return jac;
}

Eigen::VectorXd weber(const Eigen::VectorXd &x,
                      const std::list<Eigen::VectorXd> &ps) {
  Eigen::VectorXd value(ps.size());
  int i = 0;
  
  for ( auto elem : ps ) {
    value(i) = (x-elem).norm();
    i++;
  }
  return value;
}


Eigen::MatrixXd weberJ(const Eigen::VectorXd &x,
                      const std::list<Eigen::VectorXd> &ps) {
  
  Eigen::MatrixXd jac(ps.size(), x.size()); 
  int i = 0;
  Eigen::VectorXd frac = weber(x, ps);
  for ( auto elem : ps ) {
    jac.row(i++) = (x-elem)/frac(i);
  }
  return jac;
}