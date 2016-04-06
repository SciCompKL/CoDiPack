#include <codi.hpp>
#include <iostream>

template<typename Real>
void func(const Real* x, size_t l, Real* y) {
  y[0] = 0.0;
  y[1] = 1.0;
  for(size_t i = 0; i < l; ++i) {
    y[0] += x[i];
    y[1] *= x[i];
  }
}

int main(int nargs, char** args) {
  // Forward vector mode
  codi::RealForwardVec<5> x[5];
  codi::RealForwardVec<5> y[2];
  x[0] = 1.0;
  x[1] = 2.0;
  x[2] = 3.0;
  x[3] = 4.0;
  x[4] = 5.0;


  for(size_t i = 0; i < 5; ++i) {
    x[i].gradient()[i] = 1.0;
  }

  func(x, 5, y);

  double jacobi[5][2];
  for(size_t i = 0; i < 5; ++i) {
    jacobi[i][0] = y[0].getGradient()[i];
    jacobi[i][1] = y[1].getGradient()[i];
  }

  std::cout << "Forward vector mode:" << std::endl;
  std::cout << "f(1 .. 5) = (" << y[0] << ", " << y[1] << ")" << std::endl;
  for(size_t i = 0; i < 5; ++i) {
    std::cout << "df/dx_" << (i + 1) << " (1 .. 5) = (" << jacobi[i][0] << ", " << jacobi[i][1] << ")" << std::endl;
  }

  // Reverse vector mode
  codi::RealReverseVec<2> xR[5];
  codi::RealReverseVec<2> yR[2];
  xR[0] = 1.0;
  xR[1] = 2.0;
  xR[2] = 3.0;
  xR[3] = 4.0;
  xR[4] = 5.0;

  codi::RealReverseVec<2>::TapeType& tape = codi::RealReverseVec<2>::getGlobalTape();
  tape.setActive();

  for(size_t i = 0; i < 5; ++i) {
    tape.registerInput(xR[i]);
  }
  func(xR, 5, yR);
  tape.registerOutput(yR[0]);
  tape.registerOutput(yR[1]);

  tape.setPassive();

  yR[0].gradient()[0] = 1.0;
  yR[1].gradient()[1] = 1.0;
  tape.evaluate();

  double jacobiR[5][2];
  for(size_t i = 0; i < 5; ++i) {
    jacobiR[i][0] = xR[i].getGradient()[0];
    jacobiR[i][1] = xR[i].getGradient()[1];
  }

  std::cout << "Reverse vector mode:" << std::endl;
  std::cout << "f(1 .. 5) = (" << yR[0] << ", " << yR[1] << ")" << std::endl;
  for(size_t i = 0; i < 5; ++i) {
    std::cout << "df/dx_" << (i + 1) << " (1 .. 5) = (" << jacobiR[i][0] << ", " << jacobiR[i][1] << ")" << std::endl;
  }

  return 0;
}
