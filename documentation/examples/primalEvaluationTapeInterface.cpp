
#include <iostream>
#include <codi.hpp>

int main(int nargs, char** args) {
//! [Primal evaluation]
  using Real = codi::RealReversePrimalIndex;
  using Primal = typename Real::Real;
  using Identifier = typename Real::Identifier;
  using Tape = typename Real::Tape;

  Tape& tape = Real::getTape();

  Identifier x_i;
  Identifier y_i;

  // Recording
  {
    Real x = 10.0;

    tape.setActive();
    tape.registerInput(x);
    x_i = x.getIdentifier();

    Real y = 42.0 * x * x;

    tape.registerOutput(y);
    y_i = y.getIdentifier();
    tape.setPassive();
  }

  // Reverse evaluation
  tape.gradient(y_i) = 1.0;
  tape.evaluate();

  std::cout << "Gradient of dy/dx(10.0): " << tape.gradient(x_i) << std::endl;

  // Primal reevaluation and reverse evaluation
  for(int i = 0; i < 20; i += 1) {
    // Reset the gradietn of x
    tape.gradient(x_i) = 0.0;

    Primal x_v = i;
    tape.primal(x_i) = x_v;
    tape.evaluatePrimal();

    Primal y_v = tape.primal(y_i);
    tape.gradient(y_i) = 1.0;

    tape.evaluate();

    std::cout << "Value of f(" << x_v << ") = " << y_v << ", Gradient of df/dx(" << x_v << ") = " << tape.gradient(x_i) << std::endl;
  }
//! [Primal evaluation]
}
