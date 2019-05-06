#include "DownhillSimplex.hpp"


template<typename OBJECTIVEFUNCTION>
DownhillSimplex<OBJECTIVEFUNCTION>::DownhillSimplex(
      ObjectiveFunction & objectiveFunction,
      Scalar xtol,
      Scalar ftol ,
      size_t maxIterations ):
    _xtol(xtol),
    _ftol(ftol),
    _maxIterations(maxIterations),
    _objectiveFunction(objectiveFunction),
    _iterations(0) { }

template<typename OBJECTIVEFUNCTION>
DownhillSimplex<OBJECTIVEFUNCTION>::~DownhillSimplex( ) {}



template<typename OBJECTIVEFUNCTION>
int DownhillSimplex<OBJECTIVEFUNCTION>::minimize(InputType & input) {
  // Pre define
  double R{1},K{0.5},E{2},S{0.5},step{1};
  // Init simplex
  int N = input.size();
  std::list<InputType> vertices;
  vertices.push_back(input);
  int ng = -1;
  for (auto i = 0; i < N; i++) {
    InputType x = input;
    x(i) += ng*step;
    ng = -ng;
    vertices.push_back(x);
  }
  // Iterations
  InputType before = input;
  bool ternimate = false;
  while ( _iterations < _maxIterations) {
    _iterations++;
    // Sort points by function value 
    vertices.sort([&](const InputType &a, InputType &b) {
      if (_objectiveFunction.energy(a) > _objectiveFunction.energy(b)) {
        return false;
      }else{
        return true;
      } 
    });
    // Calculate average point except the worst point
    InputType xm = averageExceptWorst(vertices);
    // Do one times operation on downhill simplex
    iteration(std::ref(vertices), xm, std::ref(ternimate));

    
    // Stop condition check
    InputType now = vertices.front();
    input = now;
    if ( now == before ) {
      // most time the first vertex do not change
    }else if ( (now-before).norm() < _xtol ) {
      return 1;
    }else if (fabs( _objectiveFunction.energy(now) - _objectiveFunction.energy(before)) < _ftol) {
      return 2;
    }
    // Restart the simplex because it may converage to local region minimize
    // restart makes it overcome the local region minimize
    if ( fabs(_objectiveFunction.energy(vertices.front()) - _objectiveFunction.energy(vertices.back())) < _xtol) {
      restart(vertices);
    }
  }
  input = vertices.front();
  return -1;
}



template<typename OBJECTIVEFUNCTION>
typename DownhillSimplex<OBJECTIVEFUNCTION>::InputType DownhillSimplex<OBJECTIVEFUNCTION>::averageExceptWorst(const std::list<InputType> & vers) {
  int N = vers.size();
  InputType total = -vers.back();
  for( auto it = vers.begin(); it != vers.end(); it++) {
    total += *it;
  }
  return total/(N-1);
}

template<typename OBJECTIVEFUNCTION>
void DownhillSimplex<OBJECTIVEFUNCTION>::iteration(std::list<InputType> &vs, InputType xm, bool & terminate) {
  // Pre define
  double R{1},K{0.5},E{2},S{0.5},step{2};
  // Point x_(n+1)
  InputType xn1 = vs.back();
  // Point x_(r)
  InputType xr = xm+ R * (xm-xn1);
  Scalar Exr = _objectiveFunction.energy(xr);
  Scalar Ex1 = _objectiveFunction.energy(vs.front());
  int N = vs.size();
  auto it = vs.begin();
  for (int i = 0; i < (N-2); i++) {
    it++;
  }
  // Point x_(n)
  InputType xn = *it;
  Scalar Exn = _objectiveFunction.energy(xn);
  // reflection
  if ( Ex1 <= Exr && Exr < Exn ) { 
    vs.back() = xr;
    return ;
  }
  // Expansion
  if ( Exr < Ex1 ) { 
    InputType xe = xm + E * (xm - xn1);
    Scalar Exe = _objectiveFunction.energy(xe);
    if (Exe < Exr) {
      vs.back() = xe;
    }else {
      vs.back() = xr;
    }
    return ;
  }
  // Contraction
  if ( Exr >= Exn ) { 
    Scalar Exn1 = _objectiveFunction.energy(xn1);
    if (  Exr < Exn1) {
      InputType xoc = xm + K * ( xr - xm );
      Scalar Exoc = _objectiveFunction.energy(xoc); 
      if ( Exoc < Exr ) {
        vs.back() = xoc;
        return ;
      }else {
        shrink(vs);
        return ;
      }
    } else if ( Exn1 <= Exr) {
      InputType xic = xm - K * ( xm - xn1 );
      Scalar Exic = _objectiveFunction.energy(xic);
      if ( Exic < Exn1 ) {
        vs.back() = xic;
        return ;
      }else {
        shrink(vs);
        return ;
      }
    }
  }
  return ;
}

template<typename OBJECTIVEFUNCTION>
void DownhillSimplex<OBJECTIVEFUNCTION>::shrink(std::list<InputType> & vs) {
  double S = 0.5;
  InputType x1 = vs.front();
  for ( auto it = ++vs.begin(); it != vs.end(); it++ ) {
    InputType xi = *it;
    *it = x1 + S * ( xi - x1 );
  }
  return ;
}

// Restart the simplex
// Keep the best point as initial point,
// generate a new simplex based on this point
template<typename OBJECTIVEFUNCTION>
void DownhillSimplex<OBJECTIVEFUNCTION>::restart(std::list<InputType> &vs) {
  double step{1};
  InputType min = vs.front();
  int N = min.size();
  auto it = vs.begin();
  it++;
  int i = 0;
  int ng = -1;
  for ( ; it != vs.end(); it++) {
    InputType x = min;
    x(i) += ng * step;
    ng = -ng;
    *it = x;
    i++;
  }
}