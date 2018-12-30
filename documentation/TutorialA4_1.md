Tutorial A4.1: OpenMP reverse mode evaluation {#TutorialA4_1}
=============================================================

[Tutorial A4](@ref TutorialA4) introduces the new vector mode helper.
This helper can also be used to parallelize multiple reverse mode interpretations
with OpenMP.

The example from tutorial A4 with the vector helper is:
~~~~{.cpp}
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
~~~~
The reverse evaluation section will now be parallelized with OpenMP. For
this the vector helper will be no longer templated with a direction. Instead
it just uses a regular double for the adjoint vector.
Otherwise we will use two OpenMP threads for the evaluation:
~~~~{.cpp}
  double jacobiR[5][2];
  #pragma omp parallel num_threads(2)
  {
    int tid = omp_get_thread_num();
    codi::TapeVectorHelper<codi::RealReverse, double> vh;
    vh.gradient(yR[tid].getGradientData()) = 1.0;
    vh.evaluate();

    for(size_t i = 0; i < 5; ++i) {
      jacobiR[i][tid] = vh.getGradient(xR[i].getGradientData());
    }
  }
~~~~
It also possible to use OpenMP parallelization and multiple directions
together.

The full code for the tutorial is:
~~~~{.cpp}
#include <codi.hpp>

#include <omp.h>
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

void openMp() {
  std::cout << "codi::RealReverse( OpenMP):" << std::endl;
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

  double jacobiR[5][2];
  #pragma omp parallel num_threads(2)
  {
    int tid = omp_get_thread_num();
    codi::TapeVectorHelper<codi::RealReverse, double> vh;
    vh.gradient(yR[tid].getGradientData()) = 1.0;
    vh.evaluate();

    for(size_t i = 0; i < 5; ++i) {
      jacobiR[i][tid] = vh.getGradient(xR[i].getGradientData());
    }
  }

  std::cout << "Reverse vector mode:" << std::endl;
  std::cout << "f(1 .. 5) = (" << yR[0] << ", " << yR[1] << ")" << std::endl;
  for(size_t i = 0; i < 5; ++i) {
    std::cout << "df/dx_" << (i + 1) << " (1 .. 5) = (" << jacobiR[i][0] << ", " << jacobiR[i][1] << ")" << std::endl;
  }
}

int main(int nargs, char** args) {

  vectorHelper();

  codi::RealReverse::getGlobalTape().reset();
  openMp();

  return 0;
}
~~~~
