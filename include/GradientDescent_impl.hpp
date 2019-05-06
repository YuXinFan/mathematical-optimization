#include "GradientDescent.hpp"

template<typename OBJECTIVEFUNCTION>
GradientDescent<OBJECTIVEFUNCTION>::GradientDescent(
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
GradientDescent<OBJECTIVEFUNCTION>::~GradientDescent() { }

template<typename OBJECTIVEFUNCTION>
int GradientDescent<OBJECTIVEFUNCTION>::minimize(InputType &input) {
  InputType &x = std::ref(input);

  while (_iterations <= _maxIterations) {
    _iterations++;
    // Calculate Jacobian and value ValueType f(x)
    ValueType value = _objectiveFunction(x);
    JacobianType jac = _objectiveFunction.df(x);
    Eigen::Transpose<JacobianType> jac_trans = jac.transpose();
    // Calculate -gradient
    InputType h = -2*jac_trans * value;
    // Line search a coeff of gradient
    double alpha = lineSearch(x, h);

    InputType delta = alpha*h;
    // Apply step to x
    InputType x_new = x + delta;
    
    Scalar energy = _objectiveFunction.energy(x);
    Scalar energy_new = _objectiveFunction.energy(x_new);
    double energy_change = fabs(energy_new-energy);
    
    x = x_new;  
    // Stop condition check
    if ( delta.norm() < _xtol ) {  
      return 1;
    }
    if ( energy_change < _ftol ) { 
      return 2;
    }
    if ( fabs(energy_new) < 10*_ftol ) {  
      return 2;
    }  
  }
  return -1;
}

template<typename OBJECTIVEFUNCTION>
double GradientDescent<OBJECTIVEFUNCTION>::lineSearch(InputType & x, InputType &h) {
  // Pre define 
  double left{0}, right{10.f},mid{0}, tau{0.1};
  double norm = h.norm();
  // Stop condition
  double condition = fabs(lineSearchDF(x,h,0) * tau) ;
  // Find intervel ensure lineSearchDF=0 in it
  while (lineSearchDF(x, h, right) <= 0) {
    right *= 2;
  }

  // Iteration 100 times max
  int it = 0;
  while ( it < 100 ) {
    it++;
    mid = ( left + right )/2;
    double df = lineSearchDF(x, h, mid);
    // Shrink interval
    if (df == 0) {
      return mid;
    }else if (df > 0) {
      right = mid;
    }else {
      left = mid;
    }
    // Stop condition check
    if ( fabs(df) <= condition ) {
      break;
    }
    if ( fabs(df) <= 0.1) {
      break;
    }
  }
  return mid;
}

template<typename OBJECTIVEFUNCTION>
double GradientDescent<OBJECTIVEFUNCTION>::lineSearchDF(InputType & x, InputType &h, double alpha) {
  // Line search gradient calculation
  double phi = 2*h.transpose()* (_objectiveFunction.df(x+alpha*h).transpose()*_objectiveFunction(x+alpha*h));
  return phi;
} 
