
#include <iostream>
#include <codi.hpp>

int main(int nargs, char** args) {
//! [Reverse mode AD]
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

  // Reverse evaluation
  y.setGradient(1.0);
  tape.evaluate();

  std::cout << "Gradient of dy/dx: " << x.getGradient() << std::endl;
//! [Reverse mode AD]
}
