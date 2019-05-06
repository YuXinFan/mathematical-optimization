#ifndef CS133_HW56_GAUSSNEWTON_HPP_
#define CS133_HW56_GAUSSNEWTON_HPP_

#include <stdlib.h>


template< typename OBJECTIVEFUNCTION >
class GaussNewton {
 public:
  typedef OBJECTIVEFUNCTION ObjectiveFunction;
  typedef typename ObjectiveFunction::Scalar Scalar;
  typedef typename ObjectiveFunction::InputType InputType;
  typedef typename ObjectiveFunction::ValueType ValueType;
  typedef typename ObjectiveFunction::JacobianType JacobianType;
  
  GaussNewton(
      ObjectiveFunction & objectiveFunction,
      Scalar xtol = 0.0001,
      Scalar ftol = 0.0001,
      size_t maxIterations = 100 );
  virtual ~GaussNewton();

  //return value: -1 if optimization finished prematurely (stopped because maxiterations was reached)
  //               1 if finished because of norm(step) < xtol
  //               2 if finished because of norm(energyBefore-energyAfter) < ftol
  int minimize( InputType & input );
  double lineSearch( InputType x, InputType h);
  double lineSearchDF(InputType x, InputType h , double alpha);

private:
  ObjectiveFunction & _objectiveFunction;

  Scalar _xtol;
  Scalar _ftol;
  size_t _maxIterations;
  size_t _iterations;
};

#include "GaussNewton_impl.hpp"

#endif