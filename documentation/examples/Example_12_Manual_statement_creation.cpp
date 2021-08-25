//! [Example 12 - Manual statement creation]

#include <codi.hpp>
#include <iostream>
#include <algorithm>

using Real = codi::RealReverse;
using Tape = typename Real::Tape;

//! [Function]
// Evaluates w = (1 y y^2 ... y^(n-1)) A (1 x x^2 ... x^(n-1))^T (Standard 2D polynomial evaluation)
template<typename Real>
Real poly2D( const Real x, const Real y, const double* A, size_t n) {
    Real w = Real();

    Real curX = (Real)1.0;
    for(size_t i = 0; i < n; ++i) {

      Real curY = (Real)1.0;
      for(size_t j = 0; j < n; ++j) {

        w += A[i + j * n] * curX * curY;
        curY *= y;
      }

      curX *= x;
    }

    return w;
}
//! [Function]

template<typename Real>
Real poly2D_dx( const Real x, const Real y, const double* A, size_t n) {
    Real w = Real();

    Real curX = (Real)1.0;
    for(size_t i = 1; i < n; ++i) {

      Real curY = (Real)1.0;
      for(size_t j = 0; j < n; ++j) {
        w += (Real)i * A[i + j * n] * curX * curY;

        curY *= y;
      }

      curX *= x;
    }

    return w;
}

template<typename Real>
Real poly2D_dy( const Real x, const Real y, const double* A, size_t n) {
    Real w = Real();

    Real curX = (Real)1.0;
    for(size_t i = 0; i < n; ++i) {

      Real curY = (Real)1.0;
      for(size_t j = 1; j < n; ++j) {
        w += (Real)j * A[i + j * n] * curX * curY;

        curY *= y;
      }

      curX *= x;
    }

    return w;
}

void runExample(int mode) {
  Real u = 3.0;

  Tape& tape = Real::getTape();
  tape.setActive();
  tape.registerInput(u);

  double A[3*3] = { 1.0,    0.5,   0.25,
                    0.0,    1.0,   0.75,
                    0.25,   0.0,   1.0};
  Real x = cos(u);
  Real y = sin(u);

  Real o;

  // Step 1: Create the helper
  codi::StatementPushHelper<codi::RealReverse> ph;

  if(1 == mode) { // No special handling
    std::cout << "Running regular differentiation without statement handling." << std::endl;
    o = poly2D(x, y, A, 3);
  } else if(2 == mode) { // External function with primal function implementation
    std::cout << "Running differentiation with manual statement handling: seperate push of Jacobians." << std::endl;

    // Step 2: Compute the value with regular double values.
    double o_p = poly2D(x.getValue(), y.getValue(), A, 3);
    double jac[2];
    jac[0] = poly2D_dx(x.getValue(), y.getValue(), A, 3);
    jac[1] = poly2D_dy(x.getValue(), y.getValue(), A, 3);

    // Step 3: Push the statement on the tape
    ph.startPushStatement();
    ph.pushArgument(x, jac[0]);
    ph.pushArgument(y, jac[1]);
    ph.endPushStatement(o, o_p);
  } else if(3 == mode) { // External function with passive primal call
    std::cout << "Running differentiation with manual statement handling: Array push of Jacobians." << std::endl;

    // Step 2: Compute the value with regular double values.
    double o_p = poly2D(x.getValue(), y.getValue(), A, 3);
    double jac[2];
    jac[0] = poly2D_dx(x.getValue(), y.getValue(), A, 3);
    jac[1] = poly2D_dy(x.getValue(), y.getValue(), A, 3);

    codi::RealReverse input[2] = {x, y};

    // Step 3: Push the statement on the tape
    ph.pushStatement(o, o_p, input, jac, 2);
  } else {
    std::cerr << "Error: Unknown mode '" << mode << "'." << std::endl;
  }

  codi::RealReverse w = exp(o * o);

  tape.registerOutput(w);

  tape.setPassive();
  w.setGradient(1.0);

  tape.evaluate();

  std::cout << "Solution w: " << w << std::endl;
  std::cout << "Adjoint u: " << u.getGradient() << std::endl;

  tape.printStatistics();

  tape.reset();
}

int main(int nargs, char** args) {

  runExample(1);
  runExample(2);
  runExample(3);

  return 0;
}
//! [Example 12 - Manual statement creation]
