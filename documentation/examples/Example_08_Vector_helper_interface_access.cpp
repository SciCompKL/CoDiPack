//! [Example 8 - Vector helper interface access]

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

int main(int nargs, char** args) {

  using Real = codi::RealReverse;
  using Tape = typename Real::Tape;

  Real x[5];
  Real y[2];
  x[0] = 1.0;
  x[1] = 2.0;
  x[2] = 3.0;
  x[3] = 4.0;
  x[4] = 5.0;

  // Step 1: Perform a regular recording
  Tape& tape = Real::getTape();
  tape.setActive();

  for(size_t i = 0; i < 5; ++i) {
    tape.registerInput(x[i]);
  }

  func(x, 5, y);

  tape.registerOutput(y[0]);
  tape.registerOutput(y[1]);

  tape.setPassive();

  // Step 1: Create the helper and the get access to the vector interface
  codi::CustomAdjointVectorInterface<Real>* vh = new codi::CustomAdjointVectorHelper<Real, codi::Direction<double, 2> >();
  codi::VectorAccessInterface<Real::Real, Real::Identifier>* ai = vh->getVectorInterface();

  for(size_t dim = 0; dim < ai->getVectorSize(); ++dim) {
    ai->updateAdjoint(y[dim].getIdentifier(), dim, 1.0);   // Step 2: Set the seeding with the interface
  }

  vh->evaluate();                                          // Step 3: Call evaluate on the vector helper.

  // Step 4: Get the gradients from the interface.
  codi::Jacobian<double> jacobian(2,5);
  for(size_t i = 0; i < 5; ++i) {
    for(size_t dim = 0; dim < ai->getVectorSize(); ++dim) {
      jacobian(dim,i) = ai->getAdjoint(x[i].getIdentifier(), dim);
    }
  }

  std::cout << "Custom adjoint vector interface:" << std::endl;
  std::cout << "f(1 .. 5) = (" << y[0] << ", " << y[1] << ")" << std::endl;
  std::cout << "df/dx (1 .. 5) = \n" << jacobian << std::endl;

  tape.reset();

  delete vh;

  return 0;
}
//! [Example 8 - Vector helper interface access]
