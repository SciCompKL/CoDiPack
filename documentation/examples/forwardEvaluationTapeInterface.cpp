
#include <iostream>
#include <codi.hpp>

int main(int nargs, char** args) {
//! [Forward tape evaluation]
  using Real = codi::RealReverse;
  using Tape = typename Real::Tape;

  Tape& tape = Real::getTape();

  // Recording
  Real x = 10.0;

  tape.setActive();
  tape.registerInput(x);

  Real y = 42.0 * x * x;

  tape.registerOutput(y);
  tape.setPassive();

  // Forward evaluation
  x.setGradient(1.0);
  tape.evaluateForward();

  std::cout << "Gradient of dy/dx: " << y.getGradient() << std::endl;
//! [Forward tape evaluation]
}
