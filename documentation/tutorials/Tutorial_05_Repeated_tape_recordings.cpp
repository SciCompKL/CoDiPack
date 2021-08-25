//! [Tutorial 5 - Repeated tape recordings]

#include <codi.hpp>
#include <iostream>


//! [Function]
template<typename Real>
Real func(const Real& x, bool updateGlobal, Real& global) {
  Real t = x * x;
  if(updateGlobal) {
    global = t * x;
  }
  Real t2 = x + t;
  return global * t2;
}
//! [Function]

template<typename Real>
void run() {

  using Tape = typename Real::Tape;
  Real global = Real(0.0);

  Real x = 4.0;

  Tape& tape = Real::getTape();
  tape.setActive();

  // Step 1: Compute the gradient and update the global variable
  tape.registerInput(x);
  Real y = func(x, true, global);
  tape.registerOutput(y);

  tape.setPassive();
  y.setGradient(1.0);
  tape.evaluate();

  std::cout << "Update global:" << std::endl;
  std::cout << "f(4.0) = " << y << std::endl;
  std::cout << "df/dx(4.0) = " << x.getGradient() << std::endl;

  tape.reset();

  // Step 2: Compute the gradient but do not update the global variable.
  tape.setActive();

  tape.registerInput(x);
  y = func(x, false, global);
  tape.registerOutput(y);

  tape.setPassive();
  y.setGradient(1.0);
  tape.evaluate();

  std::cout << "No update global:" << std::endl;
  std::cout << "f(4.0) = " << y << std::endl;
  std::cout << "df/dx(4.0) = " << x.getGradient() << std::endl;

  tape.reset();

  // Step 3: Reset the identifier on the global variable.
  tape.setActive();

  tape.deactivateValue(global);

  tape.registerInput(x);
  y = func(x, false, global);
  tape.registerOutput(y);

  tape.setPassive();
  y.setGradient(1.0);
  tape.evaluate();

  std::cout << "No update global with reset:" << std::endl;
  std::cout << "f(4.0) = " << y << std::endl;
  std::cout << "df/dx(4.0) = " << x.getGradient() << std::endl;

  tape.reset();
}

int main(int nargs, char** args) {

  std::cout << "With linear index management:" << std::endl;
  run<codi::RealReverse>();
  std::cout << std::endl;
  std::cout << "With index reuse management:" << std::endl;
  run<codi::RealReverseIndex>();

  return 0;
}
//! [Tutorial 5 - Repeated tape recordings]
