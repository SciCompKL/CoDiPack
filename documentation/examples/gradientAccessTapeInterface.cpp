
#include <iostream>
#include <codi.hpp>

int main(int nargs, char** args) {
//! [Gradient Access]
  using Real = codi::RealReverse;
  using Tape = typename Real::Tape;
  using Identifier = typename Tape::Identifier;

  Tape& tape = Real::getTape();

  Identifier x_in;
  Identifier x_out;
  {
    Real x = 10.0;

    tape.setActive();
    tape.registerInput(x);
    x_in = x.getIdentifier(); // Identifier of x when it is defined as an input

    // Do some heavy computation
    x = 42.0 * x * x;

    tape.registerOutput(x);
    x_out = x.getIdentifier(); // Identifier of x when it is defined as an output

    tape.setPassive();
  }

  tape.setGradient(x_out, 1.0); // Use this interface to set the gradient of x when it was defined as an output
  tape.evaluate();

  std::cout << "Gradient of df/dx: " << tape.getGradient(x_in) << std::endl;
//! [Gradient Access]
}
