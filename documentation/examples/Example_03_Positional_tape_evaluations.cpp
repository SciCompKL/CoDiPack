//! [Example 3 - Positional tape evaluation]

#include <iostream>
#include <codi.hpp>

//! [Function]
template<typename Real>
Real func(Real a, Real b) {
  Real x = a + b;
  Real y = a - b;
  Real w = sqrt(x * x + y * y);

  return cos(w);
}
//! [Function]

//! [Positional evaluation]
using Real = codi::RealReverse;
using Gradient = typename Real::Gradient;
using Tape = typename Real::Tape;
using Position = typename Tape::Position;

Real funcInner(Real x, Real y) {
  return sqrt(x * x + y * y);
}

void positionalExample() {

  Tape& tape = Real::getTape();

  // Recording
  Real a = 10.0;
  Real b = 4.0;

  tape.setActive();
  tape.registerInput(a);
  tape.registerInput(b);

  Real u1 = a + b;
  Real u2 = a - b;

  // Now comes a really long and complicated function. Tape it, reverse it and then store the result on the tape.
  // Step 1: Store the position
  Position begin = tape.getPosition();

  // Record the function
  Real w = funcInner(u1, u2);

  // Step 2: Reverse part of the tape
  w.gradient() = 1.0;
  tape.evaluate(tape.getPosition(), begin);
  Gradient u1_d = u1.gradient();
  Gradient u2_d = u2.gradient();

  // Step 3: Cleanup the reversal
  tape.resetTo(begin);
  u1.gradient() = Gradient();
  u2.gradient() = Gradient();

  // Now store the computed gradient data
  tape.storeManual(w.value(), w.getIdentifier(), 2);
  tape.pushJacobianManual(u1_d, u1.value(), u1.getIdentifier());
  tape.pushJacobianManual(u2_d, u2.value(), u2.getIdentifier());

  Real y = cos(w);

  tape.registerOutput(y);
  tape.setPassive();

  // Reverse evaluation
  y.setGradient(1.0);
  tape.evaluate();

  std::cout << "Positional example:" << std::endl;
  std::cout << "Gradient of dy/da: " << a.getGradient() << std::endl;
  std::cout << "Gradient of dy/db: " << b.getGradient() << std::endl;

  tape.reset();
}
//! [Positional evaluation]

void validation() {

  Tape& tape = Real::getTape();

  // Recording
  Real a = 10.0;
  Real b = 4.0;

  tape.setActive();
  tape.registerInput(a);
  tape.registerInput(b);

  Real y = func(a, b);

  tape.registerOutput(y);
  tape.setPassive();

  y.setGradient(1.0);
  tape.evaluate();

  std::cout << "Validation:" << std::endl;
  std::cout << "Gradient of dy/da: " << a.getGradient() << std::endl;
  std::cout << "Gradient of dy/db: " << b.getGradient() << std::endl;

  tape.reset();
}

int main(int nargs, char** args) {
  positionalExample();
  std::cout << std::endl;
  validation();
}

//! [Example 3 - Positional tape evaluation]
