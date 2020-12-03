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

  using Real = codi::RealReverse;
  using Tape = typename Real::Tape;

  Real xR[5];
  Real yR[2];
  xR[0] = 1.0;
  xR[1] = 2.0;
  xR[2] = 3.0;
  xR[3] = 4.0;
  xR[4] = 5.0;

  Tape& tape = Real::getGlobalTape();
  tape.setActive();

//! [Custom Adjoint Vector Helper]
  for(size_t i = 0; i < 5; ++i) {
    tape.registerInput(xR[i]);
  }
  func(xR, 5, yR);
  tape.registerOutput(yR[0]);
  tape.registerOutput(yR[1]);

  tape.setPassive();

  codi::CustomAdjointVectorHelper<Real, codi::Direction<double, 2> > vh;
  vh.gradient(yR[0].getIdentifier())[0] = 1.0;
  vh.gradient(yR[1].getIdentifier())[1] = 1.0;
  vh.evaluate();

  double jacobiR[5][2];
  for(size_t i = 0; i < 5; ++i) {
    jacobiR[i][0] = vh.getGradient(xR[i].getIdentifier())[0];
    jacobiR[i][1] = vh.getGradient(xR[i].getIdentifier())[1];
  }
//! [Custom Adjoint Vector Helper]

  std::cout << "Reverse vector mode:" << std::endl;
  std::cout << "f(1 .. 5) = (" << yR[0] << ", " << yR[1] << ")" << std::endl;
  for(size_t i = 0; i < 5; ++i) {
    std::cout << "df/dx_" << (i + 1) << " (1 .. 5) = (" << jacobiR[i][0] << ", " << jacobiR[i][1] << ")" << std::endl;
  }

  return 0;
}
