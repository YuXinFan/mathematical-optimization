#ifndef CS133_HW56_GRADIENTDESCENT_HPP_
#define CS133_HW56_GRADIENTDESCENT_HPP_

#include <stdlib.h>


template< typename OBJECTIVEFUNCTION >
class GradientDescent {
 public:
  typedef OBJECTIVEFUNCTION ObjectiveFunction;
  typedef typename ObjectiveFunction::Scalar Scalar;
  typedef typename ObjectiveFunction::InputType InputType;
  typedef typename ObjectiveFunction::ValueType ValueType;
  typedef typename ObjectiveFunction::JacobianType JacobianType;
  
  GradientDescent(
      ObjectiveFunction & objectiveFunction,
      Scalar xtol = 0.00001,
      Scalar ftol = 0.00001,
      size_t maxIterations = 10000 );
  virtual ~GradientDescent();

  //return value: -1 if optimization finished prematurely (stopped because maxiterations was reached)
  //               1 if finished because of norm(step) < xtol
  //               2 if finished because of norm(energyBefore-energyAfter) < ftol
  int minimize( InputType & input );

  double lineSearch( InputType & input, InputType &grad);
  double lineSearchDF(InputType & input, InputType &grad, double alpha);
private:
  ObjectiveFunction & _objectiveFunction;

  Scalar _xtol;
  Scalar _ftol;
  size_t _maxIterations;
  size_t _iterations;
};

#include "GradientDescent_impl.hpp"

#endif