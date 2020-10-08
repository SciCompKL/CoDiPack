
#include <iostream>
#include <codi.hpp>

//! [Positional evaluation]
using Real = codi::RealReverse;
using Gradient = typename Real::Gradient;
using Tape = typename Real::Tape;
using Position = typename Tape::Position;

Real func(Real x, Real y) {
  return sqrt(x * x + y * y);
}

int main(int nargs, char** args) {

  Tape& tape = Real::getGlobalTape();

  // Recording
  Real a = 10.0;
  Real b = 4.0;

  tape.setActive();
  tape.registerInput(a);
  tape.registerInput(b);

  Real u1 = a + b;
  Real u2 = a - b;

  // Now comes a really long and complicated function. Tape it, reverse it and then store the result on the tape.
  Position begin = tape.getPosition();

  // Record the function
  Real w = func(u1, u2);

  // Reverse it
  w.gradient() = 1.0;
  tape.evaluate(tape.getPosition(), begin);
  Gradient u1_d = u1.gradient();
  Gradient u2_d = u2.gradient();

  // Cleanup
  tape.resetTo(begin);
  u1.gradient() = Gradient();
  u2.gradient() = Gradient();

  // Now store the computed gradient data
  tape.storeManual(w.value(), w.getIdentifier(), 2);
  tape.pushJacobiManual(u1_d, u1.value(), u1.getIdentifier());
  tape.pushJacobiManual(u2_d, u2.value(), u2.getIdentifier());

  Real y = cos(w);

  tape.registerOutput(y);
  tape.setPassive();

  // Reverse evaluation
  y.setGradient(1.0);
  tape.evaluate();

  std::cout << "Gradient of dy/da: " << a.getGradient() << std::endl;
  std::cout << "Gradient of dy/db: " << b.getGradient() << std::endl;
}
//! [Positional evaluation]
