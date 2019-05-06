#include "LevenbergMarquardt.hpp"

template<typename OBJECTIVEFUNCTION>
LevenbergMarquardt<OBJECTIVEFUNCTION>::LevenbergMarquardt(
    ObjectiveFunction & objectiveFunction,
    Scalar xtol ,
    Scalar ftol ,
    size_t maxIterations ):
    _xtol(xtol),
    _ftol(ftol),
    _maxIterations(maxIterations),
    _objectiveFunction(objectiveFunction),
    _iterations(0) { }

template<typename OBJECTIVEFUNCTION>
LevenbergMarquardt<OBJECTIVEFUNCTION>::~LevenbergMarquardt() { }

template<typename OBJECTIVEFUNCTION>
int LevenbergMarquardt<OBJECTIVEFUNCTION>::minimize(InputType &input) {
  // Pre define
  InputType &x= std::ref(input);
  Eigen::MatrixXd B(input.size(), input.size());
  Eigen::MatrixXd I;
  I.setIdentity(input.size(), input.size());
  double u{1.f}, tau{0.1f}, v{2};

  
  JacobianType jac = _objectiveFunction.df(x);
  // Hesseion matrix B
  B = jac.transpose() * jac;
  ValueType value = _objectiveFunction(x); 
  // Calculate init gradient
  InputType h = jac.transpose() * value;
  // Calculate init u 
  u = tau * B.diagonal().maxCoeff();
  while ( _iterations <= _maxIterations ) {
    _iterations++;
    // Calculate LM step
    InputType hlm = (B + u * I).colPivHouseholderQr().solve(-h);
    
    InputType x_new = x + hlm;
    // Calculate accept condition
    double d_F = _objectiveFunction.energy(x) - _objectiveFunction.energy(x_new);
    double d_L =  0.5 * hlm.transpose() * ( u * hlm - h );
    double q = d_F / d_L;
    // If accept new x
    if ( q > 0 ) {
      Scalar energy = _objectiveFunction.energy(x);
      // Update x
      x = x_new;
      // Update B and h
      value = _objectiveFunction(x); 
      Scalar energy_new = _objectiveFunction.energy(x);
      jac = _objectiveFunction.df(x);
      B = jac.transpose() * jac;
      h = jac.transpose() * value;
      // Stop condition check
      if ( h.norm() < 0.0001 ) {
        return 1;
      }
      if ( fabs(energy_new - energy) < _ftol ) {
        return 2;
      }
      if (hlm.norm() < _xtol ) {
        return 1;
      }
      // Update trust region size 
      double update = (1-(2*q-1)*(2*q-1)*(2*q-1));
      if ( update > (1/3.f)) {
        u = u * update;
      }else {
        u = u * 1/3.f;
      }
      v = 2; 
    // If not accept, keep x and update trust region size
    }else {
      u = u * v;
      v = v * 2;
    }
  }
  return -1;
}
