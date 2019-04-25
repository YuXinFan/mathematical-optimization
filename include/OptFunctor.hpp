#ifndef CS133_HW56_OPTFUNCTOR_HPP_
#define CS133_HW56_OPTFUNCTOR_HPP_

#include <Eigen/Eigen>
#include <iostream>


template<typename _Scalar, int NX=Eigen::Dynamic, int NY=Eigen::Dynamic>
struct OptFunctor {  
  enum
  {
    InputsAtCompileTime = NX,
    ValuesAtCompileTime = NY
  };
  
  typedef _Scalar Scalar;
  typedef Eigen::Matrix<Scalar,InputsAtCompileTime,1> InputType;
  typedef Eigen::Matrix<Scalar,ValuesAtCompileTime,1> ValueType;
  typedef Eigen::Matrix<Scalar,ValuesAtCompileTime,InputsAtCompileTime> JacobianType;

  const int m_inputs;
  const int m_values;

  OptFunctor() :
      m_inputs(InputsAtCompileTime),
      m_values(ValuesAtCompileTime) {}
  OptFunctor(int inputs, int values) :
      m_inputs(inputs),
      m_values(values) {}

  int inputs() const {
    return m_inputs;
  }
  int values() const {
    return m_values;
  }

  virtual ValueType operator()( const InputType & x ) = 0;
  virtual JacobianType df( const InputType & x ) {
    std::cout << "Error: you tried to use analytic gradient but it has not been implemented!\n";
    JacobianType res;
    return res;
  }

  Scalar energy( const InputType & x ) {
    ValueType r = this->operator()(x);
    return r.dot(r);
  }
};

#endif