//! [Example 7 - Primal tape evaluation]

#include <codi.hpp>
#include <iostream>

using Real = codi::RealReversePrimal;  // Step 1: Use a primal value taping approach
using Tape = typename Real::Tape;

//! [Function]
Real func(const Real& x) {
  return x * x * x;
}
//! [Function]

int main(int nargs, char** args) {
  Real x = 4.0;

  Tape& tape = Real::getTape();
  // Step 2: Do a normal recording and evaluation
  tape.setActive();

  tape.registerInput(x);
  Real y = func(x);
  tape.registerOutput(y);

  tape.setPassive();

  y.setGradient(1.0);
  tape.evaluate();

  std::cout << "f(4.0) = " << y << std::endl;
  std::cout << "df/dx(4.0) = " << x.getGradient() << std::endl;

  tape.clearAdjoints();

  tape.setPrimal(x.getIdentifier(), 10.0);
  tape.evaluatePrimal();

  y.setGradient(1.0);
  tape.evaluate();

  std::cout << "f(3.0) = " << tape.getPrimal(y.getIdentifier()) << std::endl;
  std::cout << "df/dx(3.0) = " << x.getGradient() << std::endl;

  tape.reset();

  return 0;
}
//! [Example 7 - Primal tape evaluation]
