#ifndef CS133_HW56_ALLFUNCTORS_HPP_
#define CS133_HW56_ALLFUNCTORS_HPP_

#include "OptFunctor.hpp"

struct RosenbrockFunctor : public OptFunctor<double,2,2> {
 public:
  typedef typename OptFunctor<double,2,2>::Scalar Scalar;
  typedef typename OptFunctor<double,2,2>::InputType InputType;
  typedef typename OptFunctor<double,2,2>::ValueType ValueType;
  typedef typename OptFunctor<double,2,2>::JacobianType JacobianType;

  RosenbrockFunctor();
  virtual ~RosenbrockFunctor();

  virtual ValueType operator()( const InputType & x );
  virtual JacobianType df( const InputType & x );
};

struct PowellFunctor : public OptFunctor<double,4,4> {
 public:
  typedef typename OptFunctor<double,4,4>::Scalar Scalar;
  typedef typename OptFunctor<double,4,4>::InputType InputType;
  typedef typename OptFunctor<double,4,4>::ValueType ValueType;
  typedef typename OptFunctor<double,4,4>::JacobianType JacobianType;

  PowellFunctor();
  virtual ~PowellFunctor();

  virtual ValueType operator()( const InputType & x );
  virtual JacobianType df( const InputType & x );
};

struct WeberFunctor : public OptFunctor<double> {
 public:
  typedef typename OptFunctor<double>::Scalar Scalar;
  typedef typename OptFunctor<double>::InputType InputType;
  typedef typename OptFunctor<double>::ValueType ValueType;
  typedef typename OptFunctor<double>::JacobianType JacobianType;

  WeberFunctor( const std::list<Eigen::VectorXd> & ps );
  virtual ~WeberFunctor();

  virtual ValueType operator()( const InputType & x );
  virtual JacobianType df( const InputType & x );

 private:
  const std::list<Eigen::VectorXd> & _ps;
};

#endif