//! [Example 14 - ReferenceActiveType]

#include <codi.hpp>
#include <iostream>
#include <sstream>

//! [Function]
template<typename Real>
void func(Real const& x, Real& y) {
  y = 3.0*x*x*x*x + 5.0*x*x*x - 3.0*x*x + 2.0*x -4.0;
}
//! [Function]

//! [FunctionRef]
template<typename Real>
void funcRef(Real& x, Real& y) {
  codi::ReferenceActiveType<Real> xRef = x;

  y = 3.0*xRef*xRef*xRef*xRef + 5.0*xRef*xRef*xRef - 3.0*xRef*xRef + 2.0*xRef -4.0;
}
//! [FunctionRef]

int main(int nargs, char** args) {

  using Real = codi::RealReverse;
  using Tape = typename Real::Tape;

  Real x = 3.14;
  Real y;

  Tape& tape = Real::getTape();

  std::cout << "Func with standard codi type." << std::endl;
  tape.setActive();

  tape.registerInput(x);
  func(x, y);
  tape.registerOutput(y);

  tape.setPassive();
  std::cout << "f(3.14) = (" << y << ")" << std::endl;

  y.setGradient(1.0);
  tape.evaluate();

  std::cout << "df/dx = (" << x.getGradient() << ")" << std::endl;

  std::ostringstream run1;
  tape.printStatistics(run1);

  tape.reset();

  std::cout << "Func with reference codi type." << std::endl;
  tape.setActive();

  tape.registerInput(x);
  funcRef(x, y);
  tape.registerOutput(y);

  tape.setPassive();
  std::cout << "f(3.14) = (" << y << ")" << std::endl;

  y.setGradient(1.0);
  tape.evaluate();

  std::cout << "df/dx = (" << x.getGradient() << ")" << std::endl;

  std::ostringstream run2;
  tape.printStatistics(run2);

  tape.reset();

  std::cout << std::endl;
  std::cout << "Statistics for the standard codi type:" << std::endl;
  std::cout << run1.str() << std::endl << std::endl;

  std::cout << "Statistics for the reference codi type:" << std::endl;
  std::cout << run2.str() << std::endl << std::endl;

  return 0;
}
//! [Example 14 - ReferenceActiveType]
