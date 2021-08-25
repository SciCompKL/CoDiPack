//! [Example 9 - OpenMP reverse evaluation]

#include <codi.hpp>
#include <iostream>
#include <omp.h>

//! [Function]
template<typename Real>
void func(const Real* x, size_t l, Real* y) {
  y[0] = 0.0;
  y[1] = 1.0;
  for(size_t i = 0; i < l; ++i) {
    y[0] += x[i];
    y[1] *= x[i];
  }
}
//! [Function]

int main(int nargs, char** args) {

  using Real = codi::RealReverse;
  using Tape = typename Real::Tape;

  Real x[5];
  Real y[2];
  x[0] = 1.0;
  x[1] = 2.0;
  x[2] = 3.0;
  x[3] = 4.0;
  x[4] = 5.0;

  // Step 1: Perform a regular recording
  Tape& tape = Real::getTape();
  tape.setActive();

  for(size_t i = 0; i < 5; ++i) {
    tape.registerInput(x[i]);
  }

  func(x, 5, y);

  tape.registerOutput(y[0]);
  tape.registerOutput(y[1]);

  tape.setPassive();

  codi::Jacobian<double> jacobian(2,5);

  #pragma omp parallel num_threads(2)
  {
    int tid = omp_get_thread_num();
    codi::CustomAdjointVectorHelper<Real, double> vh;   // Step 1: Create the vector helper for each thread

    // Step 2: Perform a regular vector helper reverse evaluation
    vh.gradient(y[tid].getIdentifier()) = 1.0;
    vh.evaluate();
    for(size_t i = 0; i < 5; ++i) {
      jacobian(tid,i) = vh.getGradient(x[i].getIdentifier());
    }
  }

  std::cout << "Custom adjoint vector helper:" << std::endl;
  std::cout << "f(1 .. 5) = (" << y[0] << ", " << y[1] << ")" << std::endl;
  std::cout << "df/dx (1 .. 5) = \n" << jacobian << std::endl;

  tape.reset();
}
//! [Example 9 - OpenMP reverse evaluation]
