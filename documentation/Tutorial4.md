Tutorial 4: Jacobi evaluation with the reverse mode      {#Tutorial4}
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
  \tag{T4.1}
  \frac{\d f_1}{\d x_j}(x) = 1.0 \eqdot
\f]
\f[
  \tag{T4.2}
  \frac{\d f_2}{\d x_j}(x) = \prod_{i \not= j} x_i \eqdot
\f]
The derivative computation with the reverse mode of CoDiPack needs now the same seven steps as in [Tutorial 2](@ref Tutorial2):
  - An implementation of the function with the CoDiPack reverse type [RealReverse](@ref codi::RealReverse).
  - A definition of the input variables.
  - Record the function evaluation on the tape.
  - A definition of the output variables.
  - The direction of the derivative has to be set.
  - The tape recorded by CoDiPack has to be evaluated.
  - The result of the derivative has to be received.

In order to get the full Jacobi matrix the last three steps have to be evaluated multiple times.

The implementation with the RealReverse type is similar as in tutorial 2.
The type in the function is exchanged and the same is done for the driver function.
The only change is the registration of the full input and output vectors on the tape.

The function looks like:
~~~~{.cpp}
    void func(const codi::RealReverse* x, size_t l, codi::RealReverse* y) {
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
   codi::RealReverse x[5];
   codi::RealReverse y[2];
   x[0] = 1.0;
   x[1] = 2.0;
   x[2] = 3.0;
   x[3] = 4.0;
   x[4] = 5.0;

   codi::RealReverse::TapeType& tape = codi::RealReverse::getGlobalTape();
   tape.setActive();

   for(size_t i = 0; i < 5; ++i) {
     tape.registerInput(x[i]);
   }
   func(x, 5, y);
   tape.registerOutput(y[0]);
   tape.registerOutput(y[1]);

   tape.setPassive();
~~~~

The tape for \f$f\f$ is now recorded.
The equation for the AD reverse mode
\f[
  \bar x = \frac{\d f}{\d x}^T(x) \cdot \bar y \quad .
\f]
just states the multiplication of the transposed Jacobi matrix with the vector \f$\bar y \in \R^m\f$.
Because \f$m=2\f$ in this example the matrix can not be computed in one sweep.
First \f$\bar y\f$ has to be set to \f$(1.0, 0.0)\f$ and then to \f$(0.0, 1.0)\f$.
Because the tape implementation in CoDiPack keeps all adjoint variables in place,
the adjoint variables have to be reset after each evaluation.
Otherwise the adjoint variables and thus the rows of the Jacobi matrix would be summed together.
The two evaluations and the reset of the adjoint variables is written as

~~~~{.cpp}
     double jacobi[5][2];

     y[0].setGradient(1.0);
     tape.evaluate();
     for(size_t i = 0; i < 5; ++i) {
       jacobi[i][0] = x[i].getGradient();
     }

     tape.clearAdjoints();

     y[1].setGradient(1.0);
     tape.evaluate();
     for(size_t i = 0; i < 5; ++i) {
       jacobi[i][1] = x[i].getGradient();
     }
~~~~

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

    void func(const codi::RealReverse* x, size_t l, codi::RealReverse* y) {
      y[0] = 0.0;
       y[1] = 1.0;
      for(size_t i = 0; i < l; ++i) {
        y[0] += x[i];
        y[1] *= x[i];
      }
    }

    int main(int nargs, char** args) {
      codi::RealReverse x[5];
      codi::RealReverse y[2];
      x[0] = 1.0;
      x[1] = 2.0;
      x[2] = 3.0;
      x[3] = 4.0;
      x[4] = 5.0;

      codi::RealReverse::TapeType& tape = codi::RealReverse::getGlobalTape();
      tape.setActive();

      for(size_t i = 0; i < 5; ++i) {
        tape.registerInput(x[i]);
      }
      func(x, 5, y);
      tape.registerOutput(y[0]);
      tape.registerOutput(y[1]);

      tape.setPassive();
      std::cout << "f(1 .. 5) = (" << y[0] << ", " << y[1] << ")" << std::endl;

      y[0].setGradient(1.0);
      tape.evaluate();

      std::cout << "df_1/dx(1 .. 5) = (";
      for(size_t i = 0; i < 5; ++i) {
        if(0 != i) {
          std::cout << ", ";
        }
        std::cout << x[i].getGradient();
      }
      std::cout << ")" << std::endl;

      tape.clearAdjoints();
      y[1].setGradient(1.0);
      tape.evaluate();

      std::cout << "df_2/dx(1 .. 5) = (";
      for(size_t i = 0; i < 5; ++i) {
        if(0 != i) {
          std::cout << ", ";
        }
        std::cout << x[i].getGradient();
      }
      std::cout << ")" << std::endl;

      return 0;
    }
~~~~
