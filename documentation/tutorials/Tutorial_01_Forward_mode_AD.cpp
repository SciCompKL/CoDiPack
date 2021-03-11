//! [Tutorial 1 - Forward mode AD]

#include <codi.hpp>
#include <iostream>

using Real = codi::RealForward;

//! [Function]
Real func(const Real& x) {
  return x * x * x;
}
//! [Function]

int main(int nargs, char** args) {
  Real x = 4.0;

  x.setGradient(1.0);      // Step 1: Set tangent seeding
  Real y = func(x);        // Step 2: Evaluate function

  // Step 3: Access gradients
  std::cout << "f(4.0) = " << y << std::endl;
  std::cout << "df/dx(4.0) = " << y.getGradient() << std::endl;

  return 0;
}
//! [Tutorial 1 - Forward mode AD]
