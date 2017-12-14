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

void vectorType() {
  std::cout << "codi::RealReverse( vector type ):" << std::endl;
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
}

void vectorHelper() {
  std::cout << "codi::RealReverse( vector helper):" << std::endl;
  // Reverse vector mode

  codi::RealReverse xR[5];
  codi::RealReverse yR[2];
  xR[0] = 1.0;
  xR[1] = 2.0;
  xR[2] = 3.0;
  xR[3] = 4.0;
  xR[4] = 5.0;

  codi::RealReverse::TapeType& tape = codi::RealReverse::getGlobalTape();
  tape.setActive();

  for(size_t i = 0; i < 5; ++i) {
    tape.registerInput(xR[i]);
  }
  func(xR, 5, yR);
  tape.registerOutput(yR[0]);
  tape.registerOutput(yR[1]);

  tape.setPassive();

  codi::TapeVectorHelper<codi::RealReverse, codi::Direction<double, 2> > vh;
  vh.gradient(yR[0].getGradientData())[0] = 1.0;
  vh.gradient(yR[1].getGradientData())[1] = 1.0;
  vh.evaluate();

  double jacobiR[5][2];
  for(size_t i = 0; i < 5; ++i) {
    jacobiR[i][0] = vh.getGradient(xR[i].getGradientData())[0];
    jacobiR[i][1] = vh.getGradient(xR[i].getGradientData())[1];
  }

  std::cout << "Reverse vector mode:" << std::endl;
  std::cout << "f(1 .. 5) = (" << yR[0] << ", " << yR[1] << ")" << std::endl;
  for(size_t i = 0; i < 5; ++i) {
    std::cout << "df/dx_" << (i + 1) << " (1 .. 5) = (" << jacobiR[i][0] << ", " << jacobiR[i][1] << ")" << std::endl;
  }
}

void vectorHelperInterface() {
  std::cout << "codi::RealReverse( vector helper interface):" << std::endl;
  // Reverse vector mode

  codi::RealReverse xR[5];
  codi::RealReverse yR[2];
  xR[0] = 1.0;
  xR[1] = 2.0;
  xR[2] = 3.0;
  xR[3] = 4.0;
  xR[4] = 5.0;

  codi::RealReverse::TapeType& tape = codi::RealReverse::getGlobalTape();
  tape.setActive();

  for(size_t i = 0; i < 5; ++i) {
    tape.registerInput(xR[i]);
  }
  func(xR, 5, yR);
  tape.registerOutput(yR[0]);
  tape.registerOutput(yR[1]);

  tape.setPassive();

  codi::TapeVectorHelperInterface<codi::RealReverse>* vh = new codi::TapeVectorHelper<codi::RealReverse, codi::Direction<double, 2> >();
  codi::AdjointInterface<codi::RealReverse::Real>* ai = vh->getAdjointInterface();

  for(size_t dim = 0; dim < ai->getVectorSize(); ++dim) {
    ai->updateAdjoint(yR[dim].getGradientData(), dim, 1.0);
  }
  vh->evaluate();

  double jacobiR[5][2];
  for(size_t i = 0; i < 5; ++i) {
    for(size_t dim = 0; dim < ai->getVectorSize(); ++dim) {
      jacobiR[i][dim] = ai->getAdjoint(xR[i].getGradientData(), dim);
    }
  }

  delete vh;

  std::cout << "Reverse vector mode:" << std::endl;
  std::cout << "f(1 .. 5) = (" << yR[0] << ", " << yR[1] << ")" << std::endl;
  for(size_t i = 0; i < 5; ++i) {
    std::cout << "df/dx_" << (i + 1) << " (1 .. 5) = (" << jacobiR[i][0] << ", " << jacobiR[i][1] << ")" << std::endl;
  }
}

int main(int nargs, char** args) {

  vectorType();

  codi::RealReverse::getGlobalTape().reset();
  vectorHelper();

  codi::RealReverse::getGlobalTape().reset();
  vectorHelperInterface();

  return 0;
}
