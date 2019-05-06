#ifndef CS133_HW56_DHSIMPLEX_HPP_
#define CS133_HW56_DHSIMPLEX_HPP_

#include <stdlib.h>


template< typename OBJECTIVEFUNCTION >
class DownhillSimplex {
 public:
  typedef OBJECTIVEFUNCTION ObjectiveFunction;
  typedef typename ObjectiveFunction::Scalar Scalar;
  typedef typename ObjectiveFunction::InputType InputType;
  typedef typename ObjectiveFunction::ValueType ValueType;
  typedef typename ObjectiveFunction::JacobianType JacobianType;
  
  DownhillSimplex(
      ObjectiveFunction & objectiveFunction,
      Scalar xtol = 0.0001,
      Scalar ftol = 0.0001,
      size_t maxIterations = 100 );
  virtual ~DownhillSimplex();

  //return value: -1 if optimization finished prematurely (stopped because maxiterations was reached)
  //               1 if finished because of norm(step) < xtol
  //               2 if finished because of norm(energyBefore-energyAfter) < ftol
  int minimize( InputType & input );

  bool sortByEnergy(InputType &a, InputType &b);
  InputType averageExceptWorst(const std::list<InputType> & vers);
  void iteration(std::list<InputType> &vers, InputType xm, bool & terminate);
  void shrink(std::list<InputType> & vs);
  void restart(std::list<InputType> & vs);

private:
  ObjectiveFunction & _objectiveFunction;

  Scalar _xtol;
  Scalar _ftol;
  size_t _maxIterations;
  size_t _iterations;
};

#include "DownhillSimplex_impl.hpp"

#endif