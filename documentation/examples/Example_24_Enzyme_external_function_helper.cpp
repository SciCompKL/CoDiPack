//! [Example 24 - Enzyme external function helper]
#include <iostream>
#include <codi.hpp>

//! [Function implementations]
template<typename Type>
Type func(Type const& v) {
  return 42.0 * v * v;
}

void func_prim(double const* x, size_t m, double* y, size_t n, codi::ExternalFunctionUserData* d) {
  y[0] = 42.0 * x[0] * x[0];
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

  // Regular computation
  y[0] = func(x);

 #if CODI_EnableEnzyme
//! [Enzyme-differentiated function]
  codi::EnzymeExternalFunctionHelper<Real> eh;

  eh.addInput(x);
  eh.addOutput(y[1]);

  eh.template callAndAddToTape<func_prim>();
//! [Enzyme-differentiated function]

//! [Enzyme-differentiated function - short]
  codi::ExternalFunctionHelper<Real> ehS;

  eh.template callAndAddToTape<func_prim>(&x, 1, &y[2], 1);
//! [Enzyme-differentiated function - short]
#else
  std::cerr << "Enzyme is not enabled for CoDiPack. Enable it with: -DCODI_EnableEnzyme=1." << std::endl;
#endif

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
//! [Example 24 - Enzyme external function helper]
