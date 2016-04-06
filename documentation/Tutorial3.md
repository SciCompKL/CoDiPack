Tutorial 3: Jacobi evaluation with the forward mode      {#Tutorial3}
============

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

The mathematical representation of the function is
\f[
  y = f(x) = \left( \sum_i x_i, \prod_i x_i \right)^T \eqdot
\f]
As the function is quite simple the Jacobi of the function can be computed by hand and is
\f[
  \tag{T3.1}
  \frac{\d f_1}{\d x_j}(x) = 1.0 \eqdot
\f]
\f[
  \tag{T3.2}
  \frac{\d f_2}{\d x_j}(x) = \prod_{i \not= j} x_i \eqdot
\f]

The derivative computation with the forward mode of CoDiPack needs now the same four steps as in [tutorial 1](@ref Tutorial1):
  - An implementation of the function with the CoDiPack forward type [RealForward](@ref codi::RealForward).
  - The direction of the derivative has to be set.
  - The evaluation of the function with the RealForward type.
  - The result of the derivative has to be received from the result.

In order to get the full Jacobi matrix the last three steps have to be evaluated multiple times.

The implementation with the RealForward type is similar as in tutorial 1.
The type in the function is exchanged and the same is done for the driver function.

The function looks like:
~~~~{.cpp}
    void func(const codi::RealForward* x, size_t l, codi::RealForward* y) {
      y[0] = 0.0;
      y[1] = 1.0;
      for(size_t i = 0; i < l; ++i) {
        y[0] += x[i];
        y[1] *= x[i];
      }
    }
~~~~

and the driver is:
~~~~{.cpp}
   codi::RealForward x[5];
   codi::RealForward y[2];
   x[0] = 1.0;
   x[1] = 2.0;
   x[2] = 3.0;
   x[3] = 4.0;
   x[4] = 5.0;

   func(x, 5, y);
~~~~

The equation for the AD forward mode
\f[
  \dot y = \frac{\d f}{\d x}(x) \cdot \dot x \quad .
\f]
just states the multiplication of the Jacobi matrix with the vector \f$\dot x \in \R^n\f$.
Because \f$n=5\f$ in this example the matrix can not be computed in one sweep.
\f$\dot x\f$ has to be set to the unit vectors in \f$\R^n\f$ e.g. \f$(1.0, 0.0, 0.0, 0.0, 0.0)\f$, \f$(0.0, 1.0, 0.0, 0.0, 0.0)\f$ and
the function needs to be evaluated for each unit vector.
The code for the evaluation is

~~~~{.cpp}
     double jacobi[5][2];

     for(size_t i = 0; i < 5; ++i) {
       x[i].setGradient(1.0);

       func(x, 5, y);

       jacobi[i][0] = y[0].getGradient();
       jacobi[i][1] = y[1].getGradient();

       x[i].setGradient(0.0);
     } .
~~~~
The function is evaluated five times and each time gradient of one input value is set to \f$1.0\f$.

The result for the primal computation will be \f$(15, 120)\f$ and the Jacobi will be

\f[
  \left( \begin{array}{ccccc}
  1.0 & 1.0 & 1.0 & 1.0 & 1.0 \\
  120 & 60 & 40 & 30 & 24\\
  \end{array} \right) \eqdot
\f]

The full code of the example is now:
~~~~{.cpp}
    #include <codi.hpp>
    #include <iostream>

    void func(const codi::RealForward* x, size_t l, codi::RealForward* y) {
      y[0] = 0.0;
       y[1] = 1.0;
      for(size_t i = 0; i < l; ++i) {
        y[0] += x[i];
        y[1] *= x[i];
      }
    }

    int main(int nargs, char** args) {
      codi::RealForward x[5];
      codi::RealForward y[2];
      x[0] = 1.0;
      x[1] = 2.0;
      x[2] = 3.0;
      x[3] = 4.0;
      x[4] = 5.0;

      for(size_t i = 0; i < 5; ++i) {
        x[i].setGradient(1.0);

        func(x, 5, y);

        if(0 == i) {
          std::cout << "f(1 .. 5) = (" << y[0] << ", " << y[1] << ")" << std::endl;
        }
        std::cout << "df/dx_" << (i + 1) << " (1 .. 5) = (" << y[0].getGradient() << ", " << y[1].getGradient() << ")" << std::endl;

        x[i].setGradient(0.0);
      }

      return 0;
    }
~~~~
