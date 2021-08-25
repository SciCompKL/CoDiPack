//! [Example 15 - Preaccumulation of code parts]

#include <algorithm>
#include <iostream>

#include <codi.hpp>

using Real = codi::RealReverse;
using Tape = typename Real::Tape;

//! [Function]
// time step ode in a explicit euler sheme
// x'(t) = Ax(t)
// x_n = x_c + dt * Ax(t)
void ode(const Real* start, Real* end, int steps, Real* A, double dt, size_t n) {
  Real* cur = new Real[n];

  for(size_t i = 0; i < n; ++i) {
    end[i] = start[i];
  }

  for(int curStep = 0; curStep < steps; ++curStep) {

    std::swap(end, cur);

    for(size_t i = 0; i < n; ++i) {
      end[i] = 0.0;
      for(size_t j = 0; j < n; ++j) {
        end[i] += A[j + i * n] * cur[j];
      }

      end[i] = cur[i] + dt * end[i];
    }
  }

  // we need to copy the result again if the number of steps is uneven
  if(steps % 2 == 1) {
    std::swap(end, cur);

    for(size_t i = 0; i < n; ++i) {
      end[i] = cur[i];
    }
  }

  delete[] cur;
}
//! [Function]

void compute(bool performPreAcc) {
  Real u = 3.0;

  Tape& tape = Real::getTape();
  tape.setActive();
  tape.registerInput(u);

//! [Preaccumulation region]
  Real A[4] = {u * 1.0, 0.5,
               0.0, u * -1.0};
  Real start[2] = {u * 10.0, u * 20.0};
  Real end[2];

  codi::PreaccumulationHelper<Real> ph;             // Step 1: Create the helper structure
  if(performPreAcc) {
    ph.start(start[0], start[1]);                   // Step 2: Start the preaccumulation region and specify inputs
    for(size_t i = 0; i < 4; ++i) {
      ph.addInput(A[i]);                            // Step 3: Add additional inputs (optional)
    }
  }

  ode(start, end, 1000, A, 1.0 / 1000.0, 2);        // Step 4: Evaluate the function/region for the accumulation in a normal way

  if(performPreAcc) {
    ph.addOutput(end[1]);                           // Step 5: Add additional outputs (optional)
    ph.finish(false, end[0]);                       // Step 6: Finish the preaccumulation region and specify outputs
  }

  Real w = sqrt(end[0] * end[0] + end[1] * end[1]);
//! [Preaccumulation region]

  tape.registerOutput(w);

  tape.setPassive();
  w.setGradient(1);

  tape.evaluate();

  std::cout << "Solution w: " << w << std::endl;
  std::cout << "Adjoint u: " << u.getGradient() << std::endl;

  tape.printStatistics();
  tape.reset();
}

int main(int nargs, char** args) {

  std::cout << "Without preaccumulation:" << std::endl;
  compute(false);
  std::cout << std::endl;

  std::cout << "With preaccumulation:" << std::endl;
  compute(true);

  return 0;
}
//! [Example 15 - Preaccumulation of code parts]
