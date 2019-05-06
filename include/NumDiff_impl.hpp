#include "NumDiff.hpp"

template<typename FUNCTOR, NumDiffMode mode>
NumDiff<FUNCTOR, mode>::NumDiff(const Functor &f, Scalar epsfcn): FUNCTOR(f) {
  _epsfcn = epsfcn;
}

template<typename FUNCTOR, NumDiffMode mode>
typename NumDiff<FUNCTOR, mode>::JacobianType NumDiff<FUNCTOR, mode>::df(const InputType &x) {
  ValueType val = this->operator()(x);
  // Call df(x, val) due to Dynamic JacobianType need 
  // both InputType and ValueType to init
  JacobianType jac = df(x, val);
  return jac;
}

template<typename FUNCTOR, NumDiffMode mode>
typename NumDiff<FUNCTOR, mode>::JacobianType NumDiff<FUNCTOR, mode>::df(const InputType &x, const ValueType &val) {
  JacobianType jac(val.size(), x.size());
  
  if ( mode == Central ){
    if ( _epsfcn == 0. ) {
    _epsfcn = 1e-9;
    }
    for (auto i = 0; i < x.size(); i++) {
      InputType x_left = x;
      InputType x_right = x;
      x_left(i) -= _epsfcn;
      x_right(i) += _epsfcn;
      ValueType value_left = this->operator()(x_left);
      ValueType value_right = this->operator()(x_right);
      jac.col(i) = (value_right-value_left)/(2*_epsfcn);
    }
  }else {
    if ( _epsfcn == 0. ) {
    _epsfcn = 1e-9;
    }
    InputType x_this = x;
    ValueType value_this = this->operator()(x);
    for (auto i = 0; i < x.size(); i++) {
      InputType x_forward = x;    
      x_forward(i) += _epsfcn; 
      ValueType value_forward = this->operator()(x_forward);
      jac.col(i) = (value_forward-value_this)/(_epsfcn);
    }
  }
  return jac;
}