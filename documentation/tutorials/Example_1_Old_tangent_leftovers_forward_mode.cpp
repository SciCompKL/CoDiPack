//! [Example 1 - Old tangent leftovers forward mode]

#include <codi.hpp>
#include <iostream>

using Real = codi::RealForward;

Real global = Real(0.0);

//! [Function]
Real func(const Real& x, bool updateGlobal) {
  if(updateGlobal) {
    global = x * x;
  }
  return x * global;
}
//! [Function]

int main(int nargs, char** args) {
  Real x = 4.0;

  // Step 1: Evaluate function and update global variable
  x.setGradient(1.0);
  Real y = func(x, true);

  std::cout << "Update global:" << std::endl;
  std::cout << "f(4.0, true) = " << y << std::endl;
  std::cout << "df/dx(4.0, true) = " << y.getGradient() << std::endl << std::endl;

  // Step 2: Evaluate function without updating the global variable
  x.setGradient(1.0);
  y = func(x, false);

  std::cout << "No update global:" << std::endl;
  std::cout << "f(4.0, false) = " << y << std::endl;
  std::cout << "df/dx(4.0, false) = " << y.getGradient() << std::endl << std::endl;

  // Step 3: Reset tangent on global variable
  global.setGradient(0.0);
  x.setGradient(1.0);
  y = func(x, false);

  std::cout << "No update global with reset:" << std::endl;
  std::cout << "f(4.0, false) = " << y << std::endl;
  std::cout << "df/dx(4.0, false) = " << y.getGradient() << std::endl << std::endl;

  return 0;
}
//! [Example 1 - Old tangent leftovers forward mode]
