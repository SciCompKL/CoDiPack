//! [Example 11 - External function user data]

#include <codi.hpp>
#include <iostream>

using Real = codi::RealReverse;
using Tape = typename Real::Tape;
using Identifier = typename Real::Identifier;
using RealBase = typename Real::Real;

//! [Function]
Real func(const Real& x) {
  return x * x * x;
}
//! [Function]


void extFunc_rev(Tape* t, void* d, codi::VectorAccessInterface<typename Tape::Real, Identifier>* va) {
  codi::ExternalFunctionUserData* data = (codi::ExternalFunctionUserData*)d;

  // Step 4: Get data in the same order as it was added.
  Identifier t_i = data->getData<int>();
  double scale = data->getData<double>();

  Real t_b = va->getAdjoint(t_i, 0);

  std::cout << " Reverse: t_b = " << t_b  << ", scale = " << scale << std::endl;
}

void extFunc_del(Tape* t, void* d) {
  codi::ExternalFunctionUserData* data = (codi::ExternalFunctionUserData*)d;

  // Step 5: Delete the data
  delete data;

  std::cout << " Reset: data is deleted." <<  std::endl;
}

int main(int nargs, char** args) {
  Real x = 4.0;

  Tape& tape = Real::getTape();
  tape.setActive();

  tape.registerInput(x);
  Real t = func(x);

  codi::ExternalFunctionUserData* data = new codi::ExternalFunctionUserData();  // Step 1: Create the data object
  data->addData(t.getIdentifier());                                  // Step 2: Add data
  data->addData(0.01);

  // Step 3: Add the external function with the data
  tape.pushExternalFunction(codi::ExternalFunction<Tape>::create(extFunc_rev, data, extFunc_del));

  Real y = func(t);
  tape.registerOutput(y);

  tape.setPassive();
  y.setGradient(1.0);
  tape.evaluate();

  std::cout << "f(f(4.0)) = " << y << std::endl;
  std::cout << "d(f ○ f)/dx(4.0) = " << x.getGradient() << std::endl;

  tape.reset();

  return 0;
}
//! [Example 11 - External function user data]
