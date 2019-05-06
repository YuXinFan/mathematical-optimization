#include "GaussNewton.hpp"


template<typename OBJECTIVEFUNCTION>
GaussNewton<OBJECTIVEFUNCTION>::GaussNewton(
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
GaussNewton<OBJECTIVEFUNCTION>::~GaussNewton() { }

template<typename OBJECTIVEFUNCTION>
int GaussNewton<OBJECTIVEFUNCTION>::minimize( InputType & input ) {
  InputType & x = std::ref(input);
  while ( _iterations <= _maxIterations ) {
    _iterations++; 
    ValueType value = _objectiveFunction(x);
    JacobianType jac = _objectiveFunction.df(x);
    // Calculate step by GaussNewton
    InputType h = (jac.transpose()*jac).colPivHouseholderQr().solve(-2*jac.transpose()*value);

    // Line search a coeff of step ensure it converage fast 
    double alpha = lineSearch(x, h);
    InputType x_new = x + alpha* h;

    Scalar energy = _objectiveFunction.energy(x);
    Scalar energy_new = _objectiveFunction.energy(x_new);
    x = x_new;
    // Stop condition check
    if (fabs(energy_new-energy) < _ftol) {
      return 2;
    }
    if (energy_new < _ftol) {
      return 2;
    }
    if ((alpha * h).norm() < _xtol) {
      return 1;
    }
  }
  return -1;
}

template<typename OBJECTIVEFUNCTION>
double GaussNewton<OBJECTIVEFUNCTION>::lineSearch( InputType x, InputType h) {
  // Pre define
  double left{0}, right{1.f}, mid{0}, tau{0.1};
  double norm = h.norm();
  double condition = fabs(lineSearchDF(x,h, 0) * tau);
  while (lineSearchDF(x,h, right) <= 0) {
    right *= 2;
  }
  while (true) {
    mid = (left+right)/2;
    double df = lineSearchDF(x, h, mid);
    // Shrink interval
    if (df < 0) {
      left = mid;
    }else if (df > 0) {
      right = mid;
    }else {
      return mid;
    }
    // Stop condition check
    if (fabs(df) <= condition) {
      break;
    }
    if (fabs(df) < 0.1) {
      break;
    }
  }
  return mid;
}

template<typename OBJECTIVEFUNCTION>
double GaussNewton<OBJECTIVEFUNCTION>::lineSearchDF( InputType x, InputType h, double alpha) {
  // Line search gradient calculation
  double phi = 2*h.transpose()* (_objectiveFunction.df(x+alpha*h).transpose()*_objectiveFunction(x+alpha*h));
  return phi;
}

