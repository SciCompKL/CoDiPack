//! [Example 6 - Forward tape evaluation]

#include <codi.hpp>
#include <iostream>

using Real = codi::RealReverse;
using Tape = typename Real::Tape;

//! [Function]
Real func(const Real& x) {
  return x * x * x;
}
//! [Function]

int main(int nargs, char** args) {
  Real x = 4.0;

  Tape& tape = Real::getTape();
  // Step 1: Do a normal recording
  tape.setActive();

  tape.registerInput(x);
  Real y = func(x);
  tape.registerOutput(y);

  tape.setPassive();

  x.setGradient(1.0);      // Step 2: Seed the input values
  tape.evaluateForward();  // Step 3: Perform forward evaluation

  // Step 4: Access gradients on the output values
  std::cout << "f(4.0) = " << y << std::endl;
  std::cout << "df/dx(4.0) = " << y.getGradient() << std::endl;

  tape.reset();

  return 0;
}
//! [Example 6 - Forward tape evaluation]
