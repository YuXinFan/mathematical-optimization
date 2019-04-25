#ifndef CS133_HW56_NUMDIFF_HPP_
#define CS133_HW56_NUMDIFF_HPP_

#include <stdlib.h>


enum NumDiffMode {
  Forward,
  Central
};

template<typename FUNCTOR, NumDiffMode mode=Forward>
class NumDiff : public FUNCTOR {
 public:
  typedef FUNCTOR Functor;
  typedef typename Functor::Scalar Scalar;
  typedef typename Functor::InputType InputType;
  typedef typename Functor::ValueType ValueType;
  typedef typename Functor::JacobianType JacobianType;

  enum {
    InputsAtCompileTime = Functor::InputsAtCompileTime,
    ValuesAtCompileTime = Functor::ValuesAtCompileTime
  };

  NumDiff( const Functor & f, Scalar epsfcn = 0. );

  JacobianType df(const InputType & x);
  JacobianType df(const InputType & x, const ValueType & val);

 private:
  Scalar _epsfcn;
};

#include "NumDiff_impl.hpp"

#endif

