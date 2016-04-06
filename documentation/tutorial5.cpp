#include <codi.hpp>
#include <iostream>
#include <sstream>

void func(const codi::RealReverse& x, codi::RealReverse& y) {
  y = 3.0*x*x*x*x + 5.0*x*x*x - 3.0*x*x + 2.0*x -4.0;
}

typedef codi::ReferenceActiveReal<codi::RealReverse> RefReal;

void funcRef(const codi::RealReverse& x, codi::RealReverse& y) {
  RefReal xRef = x;

  y = 3.0*xRef*xRef*xRef*xRef + 5.0*xRef*xRef*xRef - 3.0*xRef*xRef + 2.0*xRef -4.0;
}

int main(int nargs, char** args) {
  codi::RealReverse x = 3.14;
  codi::RealReverse y;

  codi::RealReverse::TapeType& tape = codi::RealReverse::getGlobalTape();

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
}
