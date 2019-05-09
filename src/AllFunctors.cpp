#include "AllFunctors.hpp"
#include "functions.hpp"

RosenbrockFunctor::RosenbrockFunctor() {}
RosenbrockFunctor::~RosenbrockFunctor() {}

RosenbrockFunctor::ValueType RosenbrockFunctor::operator()(const RosenbrockFunctor::InputType &x) {
  return rosenbrock(x);
}

RosenbrockFunctor::JacobianType RosenbrockFunctor::df(const RosenbrockFunctor::InputType &x) {
  return rosenbrockJ(x);
}

PowellFunctor::PowellFunctor() {}
PowellFunctor::~PowellFunctor() {}

PowellFunctor::ValueType PowellFunctor::operator()(const InputType &x) {
  return powell(x);
}

PowellFunctor::JacobianType PowellFunctor::df(const InputType &x) {
  return powellJ(x);
}

WeberFunctor::WeberFunctor(const std::list<Eigen::VectorXd> &ps):
    OptFunctor(ps.front().size(),ps.size()),
    _ps(ps) {}
WeberFunctor::~WeberFunctor() {}

WeberFunctor::ValueType WeberFunctor::operator()(const InputType &x) {
  return weber(x, _ps);
}

WeberFunctor::JacobianType WeberFunctor::df(const InputType &x) {
  return weberJ(x, _ps);
}