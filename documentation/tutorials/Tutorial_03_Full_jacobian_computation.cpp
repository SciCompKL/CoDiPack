//! [Tutorial 3 - Full Jacobian computation]

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

void forwardModeJacobianComputation() {

  using Real = codi::RealForward;

  Real x[5];
  Real y[2];
  x[0] = 1.0;
  x[1] = 2.0;
  x[2] = 3.0;
  x[3] = 4.0;
  x[4] = 5.0;

  codi::Jacobian<double> jacobian(2,5);

  // Step 1: Iterate over the input dimension
  for(size_t i = 0; i < 5; ++i) {
    x[i].gradient() = 1.0; // Step 2: Set the seeding for the i-th input variable

    func(x, 5, y);  // Step 3: Evaluate the function

    // Step 3: Get the gradients from the outputs.
    for(size_t j = 0; j < 2; ++j) {
      jacobian(j,i) = y[j].getGradient();
    }

    x[i].gradient() = 0.0; // Step 4: Reset the seeding for the i-th input variable
  }

  std::cout << "Forward mode Jacobian:" << std::endl;
  std::cout << "f(1 .. 5) = (" << y[0] << ", " << y[1] << ")" << std::endl;
  std::cout << "df/dx (1 .. 5) = \n" << jacobian << std::endl;
}

void reverseModeJacobianComputation() {

  using Real = codi::RealReverse;
  using Tape = typename Real::Tape;

  Real x[5];
  Real y[2];
  x[0] = 1.0;
  x[1] = 2.0;
  x[2] = 3.0;
  x[3] = 4.0;
  x[4] = 5.0;

  codi::Jacobian<double> jacobian(2,5);

  Tape& tape = Real::getTape();
  tape.setActive();

  // Step 1: Record the tape
  for(size_t i = 0; i < 5; ++i) {
    tape.registerInput(x[i]);
  }

  func(x, 5, y);

  tape.registerOutput(y[0]);
  tape.registerOutput(y[1]);

  tape.setPassive();

  // Step 2: Iterate over the output dimension
  for(size_t i = 0; i < 2; ++i) {
    y[i].gradient() = 1.0; // Step 3: Set the seeding for the i-th output variable

    tape.evaluate();

    // Step 4: Get the gradients from the inputs.
    for(size_t j = 0; j < 5; ++j) {
      jacobian(i, j) = x[j].getGradient();
    }

    tape.clearAdjoints(); // Step 5: Clear the adjoints
  }

  std::cout << "Reverse mode Jacobian:" << std::endl;
  std::cout << "f(1 .. 5) = (" << y[0] << ", " << y[1] << ")" << std::endl;
  std::cout << "df/dx (1 .. 5) = \n" << jacobian << std::endl;

  tape.reset(false);
}

int main(int nargs, char** args) {
  forwardModeJacobianComputation();

  std::cout << std::endl;

  reverseModeJacobianComputation();

  return 0;
}

//! [Tutorial 3 - Full Jacobian computation]
