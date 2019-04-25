#ifndef CS133_HW56_FUNCTIONS_HPP_
#define CS133_HW56_FUNCTIONS_HPP_

#include <Eigen/Eigen>
#include <list>


/// Rosenbrock function:
//     f(x) = 100(y-x^2)^2 + (1-x)^2
//    f1(x) = 10(y-x^2)
//    f2(x) = 1-x
Eigen::Vector2d
rosenbrock(const Eigen::Vector2d &x);

Eigen::Matrix2d
rosenbrockJ(const Eigen::Vector2d &x);


/// Powell's function:
//     f(x) = (x_1+10*x_2)^2 +
//             5*(x_3-x_4)^2 +
//             (x_2-2*x_3)^4 +
//            10*(x_1-x_4)^4
//    f1(x) = x_1+10*x_2
//    f2(x) = sqrt(5)*(x_3-x_4)
//    f3(x) = (x_2-2*x_3)^2
//    f4(x) = sqrt(10)*(x_1-x_4)^2
Eigen::Vector4d
powell(const Eigen::Vector4d &x);

Eigen::Matrix4d
powellJ(const Eigen::Vector4d &x);


// Weber's function:
//     f(x) = sum_i norm(x-p_i)^2
Eigen::VectorXd
weber(const Eigen::VectorXd &x, const std::list<Eigen::VectorXd> &ps);

Eigen::MatrixXd
weberJ(const Eigen::VectorXd &x, const std::list<Eigen::VectorXd> &ps);

#endif
