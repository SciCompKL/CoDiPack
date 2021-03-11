//! [Example 2 - Custom adjoint vector helper]

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

  Tape& tape = Real::getGlobalTape();
  tape.setActive();

//! [Custom Adjoint Vector Helper]
  // Step 1: Record a regular tape
  for(size_t i = 0; i < 5; ++i) {
    tape.registerInput(x[i]);
  }

  func(x, 5, y);

  tape.registerOutput(y[0]);
  tape.registerOutput(y[1]);

  tape.setPassive();

  // Step 2: Use for the seeding, the evaluation and the retrieval of the gradients the vector helper.
  codi::CustomAdjointVectorHelper<Real, codi::Direction<double, 2> > vh;
  vh.gradient(y[0].getIdentifier())[0] = 1.0;
  vh.gradient(y[1].getIdentifier())[1] = 1.0;
  vh.evaluate();

  codi::Jacobian<double> jacobian(2,5);
  for(size_t i = 0; i < 5; ++i) {
    jacobian(0,i) = vh.getGradient(x[i].getIdentifier())[0];
    jacobian(1,i) = vh.getGradient(x[i].getIdentifier())[1];
  }
//! [Custom Adjoint Vector Helper]

  std::cout << "Reverse vector mode:" << std::endl;
  std::cout << "f(1 .. 5) = (" << y[0] << ", " << y[1] << ")" << std::endl;
  std::cout << "df/dx (1 .. 5) = \n" << jacobian << std::endl;

  tape.reset();

  return 0;
}
//! [Example 2 - Custom adjoint vector helper]
