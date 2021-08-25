
#include <iostream>
#include <codi.hpp>

//! [Function implementations]
template<typename Type>
Type func(Type const& v) {
  return 42.0 * v * v;
}

template<typename Type>
void func_wrap(Type const& v, Type& w) {
  w = func(v);
}

void func_prim(double const* x, size_t m, double* y, size_t n, codi::ExternalFunctionUserData* d) {
  y[0] = 42.0 * x[0] * x[0];
}

void func_rev(double const* x, double* x_b, size_t m, double const* y, double const* y_b, size_t n, codi::ExternalFunctionUserData* d) {
  x_b[0] = y_b[0] * x[0] * 84.0;
}
//! [Function implementations]

int main(int nargs, char** args) {
  using Real = codi::RealReverse;
  using Tape = typename Real::Tape;

  Tape& tape = Real::getTape();

  Real x = 10.0;
  Real y[3];

  tape.setActive();
  tape.registerInput(x);

  // Regular compuation
  y[0] = func(x);

//! [Mode 1: Implemented primal function]
  codi::ExternalFunctionHelper<Real> eh;

  eh.addInput(x);
  eh.addOutput(y[1]);

  eh.callPrimalFunc(func_prim);
  eh.addToTape(func_rev);
//! [Mode 1: Implemented primal function]

//! [Mode 2: Passive primal function]
  codi::ExternalFunctionHelper<Real> ehP(true);

  ehP.addInput(x);
  ehP.addOutput(y[2]);

  ehP.callPrimalFuncWithADType(func_wrap<Real>, x, y[2]);
  ehP.addToTape(func_rev);
//! [Mode 2: Passive primal function]

  for(int i = 0; i < 3; i += 1) {
    tape.registerOutput(y[i]);
  }
  tape.setPassive();

  for(int i = 0; i < 3; i += 1) {
    tape.clearAdjoints();
    y[i].setGradient(1.0);
    tape.evaluate();
    std::cout << "Gradient of dy[" << i << "]/dx: " << x.getGradient() << std::endl;
  }
}
