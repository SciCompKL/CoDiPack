#include <codi.hpp>
#include <iostream>

codi::RealForward func(const codi::RealForward& x) {
  return x * x * x;
}

int main(int nargs, char** args) {
  codi::RealForward x = 4.0;
  x.setGradient(1.0);

  codi::RealForward y = func(x);

  std::cout << "f(4.0) = " << y << std::endl;
  std::cout << "df/dx(4.0) = " << y.getGradient() << std::endl;

  return 0;
}
