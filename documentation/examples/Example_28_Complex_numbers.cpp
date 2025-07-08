//! [Example 28 - Complex numbers]
#include <iostream>
#include <codi.hpp>

//! [Function implementations]
template<typename Type>
Type func(Type const& v) {
  Type t = v + v.real();
  return 2 * v * t;
}
//! [Function implementations]

int main(int nargs, char** args) {
  using Real = codi::RealReverse;
  using Tape = typename Real::Tape;

  Tape& tape = Real::getTape();

  std::complex<Real> x = 10.0;
  Real w = 5.0;

  tape.setActive();

  codi::RealTraits::registerInput(w);
  codi::RealTraits::registerInput(x); // Use general registration function for complex numbers.

  // Use complex numbers as usual.
  std::complex<Real> y = func(x);
  y *= w;

  codi::RealTraits::registerOutput(y);  // Use general registration function for complex numbers.

  tape.setPassive();

  // Cast to Real* is still possible.
  Real* x_p = reinterpret_cast<Real*>(&x);
  Real* y_p = reinterpret_cast<Real*>(&y);
  for(int i = 0; i < 2; i += 1) {
    tape.clearAdjoints();
    y_p[i].setGradient(1.0);
    tape.evaluate();
    std::cout << "Gradient of dy[" << i << "]/dw: " << w.getGradient() << std::endl;
    std::cout << "Gradient of dy[" << i << "]/dx[0]: " << x_p[0].getGradient() << std::endl;
    std::cout << "Gradient of dy[" << i << "]/dx[1]: " << x_p[1].getGradient() << std::endl;
  }
}
//! [Example 28 - Complex numbers]
