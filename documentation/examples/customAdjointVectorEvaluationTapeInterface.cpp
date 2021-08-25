
#include <iostream>
#include <codi.hpp>

int main(int nargs, char** args) {
//! [Custom vector]
  using Real = codi::RealReverse;
  using Tape = typename Real::Tape;

  Tape& tape = Real::getTape();

  // Recording
  Real x = 10.0;

  tape.setActive();
  tape.registerInput(x);

  Real y1 = 42.0 * x * x;
  Real y2 = 20.0 * x * x * x;

  tape.registerOutput(y1);
  tape.registerOutput(y2);

  tape.setPassive();

  // Reverse evaluation
  size_t adjointSize = tape.getParameter(codi::TapeParameters::LargestIdentifier);
  codi::Direction<double, 2>* adjoints = new codi::Direction<double, 2>[adjointSize + 1];

  adjoints[y1.getIdentifier()] = {1.0, 0.0};
  adjoints[y2.getIdentifier()] = {0.0, 1.0};

  tape.evaluate(tape.getPosition(), tape.getZeroPosition(), adjoints); // Full tape evaluation

  std::cout << "Gradient of dy1/dx: " << adjoints[x.getIdentifier()][0] << std::endl;
  std::cout << "Gradient of dy2/dx: " << adjoints[x.getIdentifier()][1] << std::endl;

  delete [] adjoints;
//! [Custom vector]
}
