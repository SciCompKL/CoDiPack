//! [Example 10 - External function helper]

#include <codi.hpp>
#include <iostream>

using Real = codi::RealReverse;
using Tape = typename Real::Tape;

using BaseReal = typename Real::Real;

//! [Function]
template<typename Number>
void solve2(const Number* A, const Number* b, Number* x) {

  // A = a[0] a[1]  A^-1 = 1/det *  a[3] -a[1]
  //     a[2] a[3]                 -a[2]  a[0]
  Number det = A[0] * A[3] - A[1] * A[2];

  x[0] = (A[3] * b[0] - A[1] * b[1]) / det;
  x[1] = (-A[2] * b[0] + A[0] * b[1]) / det;
}
//! [Function]

void solve2_primal(const BaseReal* x, size_t m, BaseReal* y, size_t n, codi::ExternalFunctionUserData* d) {
  solve2(&x[0], &x[4], y);
}

void solve2_rev(const BaseReal* x, BaseReal* x_b, size_t m, const BaseReal* y, const BaseReal* y_b, size_t n, codi::ExternalFunctionUserData* d) {
  BaseReal ATrans[4] = {x[0], x[2], x[1], x[3]};

  BaseReal s[2];
  solve2(ATrans, y_b, s);

  // Adjoint of A (\bar A = -s*x^T) (In local terms x[0-3] = -s*y^T)
  x_b[0] = -s[0] * y[0];
  x_b[1] = -s[0] * y[1];
  x_b[2] = -s[1] * y[0];
  x_b[3] = -s[1] * y[1];

  // Adjoint of b (\bar b = s) (In local terms x[4-5] = s)
  x_b[4] = s[0];
  x_b[5] = s[1];
}

void runExample(int mode) {
  Real u = 3.0;

  Tape& tape = Real::getTape();
  tape.setActive();
  tape.registerInput(u);

  Real A[4] = { u * 1.0, 0.5,  0.25, u * -1.0};
  Real b[2] = {u * 10.0, u * 20.0};

  Real x[2];

  if(1 == mode) { // No special handling
    std::cout << "Running regular differentiation without external functions." << std::endl;
    solve2(A, b, x);
  } else if(2 == mode) { // External function with primal function implementation
    std::cout << "Running differentiation with external function, primal is called via a special function implementation." << std::endl;
    codi::ExternalFunctionHelper<codi::RealReverse> eh; // Step 1: Create the helper
    for(int i = 0; i < 4; ++i) {
      eh.addInput(A[i]);                                // Step 2: Add inputs of the function
    }
    for(int i = 0; i < 2; ++i) {
      eh.addInput(b[i]);                                // Step 2: Add inputs of the function
    }

    for(int i = 0; i < 2; ++i) {
      eh.addOutput(x[i]);                               // Step 3: Add outputs of the function
    }

    eh.callPrimalFunc(solve2_primal);                   // Step 4: Call the primal with a special implementation.
    eh.addToTape(solve2_rev);                           // Step 5: Added specialized reverse function to the tape.
  } else if(3 == mode) { // External function with passive primal call
    std::cout << "Running differentiation with external function, primal is called via a passive AD evaluation." << std::endl;
    codi::ExternalFunctionHelper<codi::RealReverse> eh(true); // Step 1: Create the helper
    for(int i = 0; i < 4; ++i) {
      eh.addInput(A[i]);                                // Step 2: Add inputs of the function
    }
    for(int i = 0; i < 2; ++i) {
      eh.addInput(b[i]);                                // Step 2: Add inputs of the function
    }

    for(int i = 0; i < 2; ++i) {
      eh.addOutput(x[i]);                               // Step 3: Add outputs of the function
    }

    eh.callPrimalFuncWithADType(solve2<Real>, A, b, x); // Step 4: Call the primal with a regular function call that is not recorded
    eh.addToTape(solve2_rev);                           // Step 5: Added specialized reverse function to the tape.
  } else {
    std::cerr << "Error: Unknown mode '" << mode << "'." << std::endl;
  }

  Real w = sqrt(x[0] * x[0] + x[1] * x[1]);

  tape.registerOutput(w);

  tape.setPassive();
  w.setGradient(1.0);

  tape.evaluate();

  std::cout << "Solution w: " << w << std::endl;
  std::cout << "Adjoint u: " << u.getGradient() << std::endl;

  tape.reset();
}

int main(int nargs, char** args) {

  runExample(1);
  runExample(2);
  runExample(3);

  return 0;
}
//! [Example 10 - External function helper]
