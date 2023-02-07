
#include <iostream>
#include <codi.hpp>

//! [Manual statement push]
using Real = codi::RealReverse;
using Tape = typename Real::Tape;

using Primal = typename Real::Real;

Primal func(Primal x, Primal y) {
  return sqrt(x * x + y * y);
}

Primal func_dx(Primal x, Primal y) {
  return x / func(x, y);
}

Primal func_dy(Primal x, Primal y) {
  return y / func(x, y);
}

int main(int nargs, char** args) {

  Tape& tape = Real::getTape();

  // Recording
  Real u1 = 10.0;
  Real u2 = 4.0;

  tape.setActive();
  tape.registerInput(u1);
  tape.registerInput(u2);

  Real w = func(u1.value(), u2.value());
  tape.storeManual(w.value(), w.getIdentifier(), 2);
  tape.pushJacobianManual(func_dx(u1.value(), u2.value()), u1.value(), u1.getIdentifier());
  tape.pushJacobianManual(func_dy(u1.value(), u2.value()), u2.value(), u2.getIdentifier());

  tape.registerOutput(w);
  tape.setPassive();

  // Reverse evaluation
  w.setGradient(1.0);
  tape.evaluate();

  std::cout << "Gradient of dy/du1: " << u1.getGradient() << std::endl;
  std::cout << "Gradient of dy/du2: " << u2.getGradient() << std::endl;
}
//! [Manual statement push]
