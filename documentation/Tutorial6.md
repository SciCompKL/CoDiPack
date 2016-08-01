Tutorial 6: Vector mode demonstration      {#Tutorial6}
============

This is the same tutorial as [Tutorial 3](@ref Tutorial3) and [Tutorial 4](@ref Tutorial4) but uses the vector mode to compute the full Jacobian matrix.
The vector mode of AD can be used to evaluate multiple tangent or adjoint directions in one sweep.

The function we want to differentiate is implemented as
~~~~{.cpp}
    void func(const double* x, size_t l, double* y) {
      y[0] = 0.0;
      y[1] = 1.0;
      for(size_t i = 0; i < l; ++i) {
        y[0] += x[i];
        y[1] *= x[i];
      }
    }
~~~~

The derivative computation with the vector mode of CoDiPack needs the same steps as in the original tutorials but instead of setting the tangent or adjoint direction on one input or output variable, the tangent and adjoint directions are set for multiple input or output variables.
Each value is set to a different unit vector e.g. \f$(1, 0, \ldots, 0)\f$, \f$(0, 1, 0, \ldots, 0)\f$, etc.

This time the function is not implemented for a specific type but templated in order to allow the evaluation with the codi::RealForwardVec type and codi::RealReverseVec type.
The function looks like:
~~~~{.cpp}
    template<typename Real>
    void func(const Real* x, size_t l, Real* y) {
      y[0] = 0.0;
      y[1] = 1.0;
      for(size_t i = 0; i < l; ++i) {
        y[0] += x[i];
        y[1] *= x[i];
      }
    }
~~~~

The equation for the AD forward mode
\f[
  \dot y = \frac{\d f}{\d x}(x) \cdot \dot x
\f]
can be extended to the AD forward vector mode
\f[
  \dot Y = \frac{\d f}{\d x}(x) \cdot \dot X
\f]
were \f$\dot X \in \R^{(n \times k)}\f$ is now a matrix.
Each column of the matrix can be set to a unit vector and then the full Jacobian of \f$f\f$ can be evaluated in one sweep.
The code for the driver is:

~~~~{.cpp}
   codi::RealForwardVec<5> x[5];
   codi::RealForwardVec<5> y[2];
   x[0] = 1.0;
   x[1] = 2.0;
   x[2] = 3.0;
   x[3] = 4.0;
   x[4] = 5.0;

   double jacobi[5][2];

   for(size_t i = 0; i < 5; ++i) {
     x[i].gradient()[i] = 1.0;
   }

   func(x, 5, y);

   for(size_t i = 0; i < 5; ++i) {
     jacobi[i][0] = y[0].getGradient()[i];
     jacobi[i][1] = y[1].getGradient()[i];
   }
~~~~

The number of input variables in this example is 5 and therefore the size for the tangent direction is also chosen as 5.
This allows to compute the Jacobian in one sweep.
If the size for the tangent direction is smaller e.g. 3 then 2 sweeps would be needed.
The gradient values are now accessed via the direct accessor `.gradient()`.
This function returns a reference to the tangent direction.
The entries of the direction can be accesses via the array operator `[]`.
In order to set the unit vectors, the i-th entry of the i-th argument is set to 1.0.

The result for the primal computation will be \f$(15, 120)\f$ and the Jacobi will be

\f[
  \left( \begin{array}{ccccc}
  1.0 & 1.0 & 1.0 & 1.0 & 1.0 \\
  120 & 60 & 40 & 30 & 24\\
  \end{array} \right) \eqdot
\f]

The same principle can be used for the reverse mode of AD.
This time the dimension for the vectors is set to 2 because there are only two output arguments.

The driver for the reverse mode is:

~~~~{.cpp}
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
~~~~

The full code of the example is:
~~~~{.cpp}
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
~~~~
