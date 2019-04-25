#ifndef CS133_HW56_LEVENBERGMARQUARDT_HPP_
#define CS133_HW56_LEVENBERGMARQUARDT_HPP_

#include <stdlib.h>


template< typename OBJECTIVEFUNCTION >
class LevenbergMarquardt {
 public:
  typedef OBJECTIVEFUNCTION ObjectiveFunction;
  typedef typename ObjectiveFunction::Scalar Scalar;
  typedef typename ObjectiveFunction::InputType InputType;
  typedef typename ObjectiveFunction::ValueType ValueType;
  typedef typename ObjectiveFunction::JacobianType JacobianType;
  
  LevenbergMarquardt(
      ObjectiveFunction & objectiveFunction,
      Scalar xtol = 0.0001,
      Scalar ftol = 0.0001,
      size_t maxIterations = 100 );
  virtual ~LevenbergMarquardt();

  //return value: -1 if optimization finished prematurely (stopped because maxiterations was reached)
  //               1 if finished because of norm(step) < xtol
  //               2 if finished because of norm(energyBefore-energyAfter) < ftol
  int minimize( InputType & input );

private:
  ObjectiveFunction & _objectiveFunction;

  Scalar _xtol;
  Scalar _ftol;
  size_t _maxIterations;
  size_t _iterations;
};

#include "LevenbergMarquardt_impl.hpp"

#endif