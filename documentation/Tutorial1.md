Tutorial 1: Forward mode      {#Tutorial1}
============


The first function we want to differentiation is implemented as
~~~~{.cpp}
    double func(const double& x) {
      return x * x * x;
    }
~~~~

The mathematical representation of the function is
\f[
  y = f(x) = x * x * x = x^3 \eqdot
\f]
As the function is quite simple the jacobi of the function can be computed by hand and is
\f[
  \tag{T1.1}
  \frac{\d f}{\d x}(x) = 3 * x * x = 3 * x^2 \eqdot
\f]
The derivative computation with the forward mode of CoDiPack needs now four things:
  - An implementation of the function with the CoDiPack forward type [RealForward](@ref codi::RealForward).
  - The direction of the derivative has to be set.
  - The evaluation of the function with the RealForward type.
  - The result of the derivative has to be received from the result.

The first thing is the implementation of the function with the forward type of CoDiPack.
This is done by replacing every occurrence of the type double in the function with the type [RealForward](@ref codi::RealForward).
The overloading implementation looks like
~~~~{.cpp}
    codi::RealForward func(const codi::RealForward& x) {
      return x * x * x;
    }
~~~~

The driver for the primal evaluation of the function is inside the main routine:
~~~~{.cpp}
    int main(int nargs, char** args) {
      double x = 4.0;

      double y = func(x);

      std::cout << "f(4.0) = " << y << std::endl;

      return 0;
    }
~~~~

The overloaded function can only be called if we exchange the double types with the RealForward type.
~~~~{.cpp}
    int main(int nargs, char** args) {
      codi::RealForward x = 4.0;

      codi::RealForward y = func(x);

      std::cout << "f(4.0) = " << y << std::endl;

      return 0;
    }
~~~~

Now the implementation of the function with the RealForward type is completed.
The next step will set the direction of the derivative.
In order to get a better feeling what this means, the forward AD formulation is repeated here.
For the function \f$y = f(x)\f$  with the input variables \f$x \in \R^n\f$ and the output variables \f$y \in \R^m\f$.
The forward mode of AD computes the equation
\f[
  \dot y = \frac{\d f}{\d x}(x) \cdot \dot x \quad .
\f]
In our example the dimensions of the vector spaces \f$\R^n\f$ and \f$\R^m\f$ is one.
The formulation contains the partial derivative of \f$f\f$ with respect to the input variable \f$x\f$,
which was derived in equation (T1.1).
AD will multiply this partial derivative with the value of \f$\dot x\f$ and the result of the operation will be stored in \f$\dot y\f$.
By setting \f$\dot x = 1\f$ the result will be the gradient of \f$f\f$ at the position \f$x\f$.
The RealForward type stores the primal value and the derivative value which is denoted by the dot above the name of the value.
The primal value can be accessed with the functions getValue(), setValue() and value().
For the gradient the same functions getGradient(), setGradient() and gradient() exist.
The value of \f$\dot x\f$ can be set with the call:
~~~~{.cpp}
      x.setGradient(1.0);
~~~~
The value of \f$\dot y\f$ can be retrieved with the call:
~~~~{.cpp}
      double y_dot = y.getGradient();
~~~~

The implementation of the main method is now extended by these two calls:
~~~~{.cpp}
    int main(int nargs, char** args) {
      codi::RealForward x = 4.0;
      x.setGradient(1.0);

      codi::RealForward y = func(x);

      std::cout << "f(4.0) = " << y << std::endl;
      std::cout << "df/dx(4.0) = " << y.getGradient() << std::endl;

      return 0;
    }
~~~~

The result for the primal computation will be 64 and the derivative will be 48.

The full code of the example is now:
~~~~{.cpp}
    #include <codi.hpp>
    #include <iostream>

    codi::RealForward func(const codi::RealForward& x) {
      return x * x * x;
    }

    int main(int nargs, char** args) {
      codi::RealForward x = 4.0;
      x.setGradient(1.0);

      codi::RealForward y = func(x);

      std::cout << "f(4.0) = " << y << std::endl;
      std::cout << "df/dx(4.0) = " << y.getGradient() << std::endl;

      return 0;
    }
~~~~
