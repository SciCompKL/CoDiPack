//! [Tutorial 4 - Vector mode AD]

#include <codi.hpp>
#include <iostream>

//! [Function]
template<typename Real>
void func(const Real* x, size_t l, Real* y) {
  y[0] = 0.0;
  y[1] = 1.0;
  for(size_t i = 0; i < l; ++i) {
    y[0] += x[i];
    y[1] *= x[i];
  }
}

//! [Function]

void forwardVectorMode() {

  using Real = codi::RealForwardVec<5>;   // Step 1: Use the vector mode type

  Real x[5];
  Real y[2];
  x[0] = 1.0;
  x[1] = 2.0;
  x[2] = 3.0;
  x[3] = 4.0;
  x[4] = 5.0;

  // Step 2: Set the seeding for each vector direction
  for(size_t i = 0; i < 5; ++i) {
    x[i].gradient()[i] = 1.0;
  }

  func(x, 5, y);

  // Step 3: Get the gradients from the outputs.
  codi::Jacobian<double> jacobian(2,5);
  for(size_t i = 0; i < 5; ++i) {
    jacobian(0,i) = y[0].getGradient()[i];
    jacobian(1,i) = y[1].getGradient()[i];
  }

  std::cout << "Forward vector mode:" << std::endl;
  std::cout << "f(1 .. 5) = (" << y[0] << ", " << y[1] << ")" << std::endl;
  std::cout << "df/dx (1 .. 5) = \n" << jacobian << std::endl;
}

void reverseVectorMode() {

  using Real = codi::RealReverseVec<2>;   // Step 1: Use the vector mode type
  using Tape = typename Real::Tape;

  Real x[5];
  Real y[2];
  x[0] = 1.0;
  x[1] = 2.0;
  x[2] = 3.0;
  x[3] = 4.0;
  x[4] = 5.0;

  Tape& tape = Real::getTape();
  tape.setActive();

  for(size_t i = 0; i < 5; ++i) {
    tape.registerInput(x[i]);
  }

  func(x, 5, y);

  tape.registerOutput(y[0]);
  tape.registerOutput(y[1]);

  tape.setPassive();

  // Step 2: Set the seeding for each vector direction
  y[0].gradient()[0] = 1.0;
  y[1].gradient()[1] = 1.0;

  tape.evaluate();

  // Step 3: Get the gradients from the inputs.
  codi::Jacobian<double> jacobian(2,5);
  for(size_t i = 0; i < 5; ++i) {
    jacobian(0,i) = x[i].getGradient()[0];
    jacobian(1,i) = x[i].getGradient()[1];
  }

  std::cout << "Reverse vector mode:" << std::endl;
  std::cout << "f(1 .. 5) = (" << y[0] << ", " << y[1] << ")" << std::endl;
  std::cout << "df/dx (1 .. 5) = \n" << jacobian << std::endl;

  tape.reset();
}

int main(int nargs, char** args) {
  forwardVectorMode();

  std::cout << std::endl;

  reverseVectorMode();

  return 0;
}

//! [Tutorial 4 - Vector mode AD]
