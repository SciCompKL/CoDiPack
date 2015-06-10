Tutorial 1: Forward mode      {#Tutorial1}
============


The first function we want to differentation in the tutorial is implemented as
~~~~{.cpp}
    double func(const double& x) {
      return x * x * x;
    }
~~~~

The mathematical representation of the function is
\f[
  y = f(x) = x * x * x \quad .
\f]
We can now easyly compute the derivative of the function by hand which yields
\f[
  \frac{\partial f}{\partial x}(x) = 3 * x * x \quad .
\f]
The derivative computation with CoDiPack needs now three things:
  - An implementation of the function with the CoDiPack forward type RealForward
  - The direction of the derivate has to be set.
  - The result of the derivative has to be recevied from the result.

The first thing is the implementation of the function with the forward type of CoDiPack.
This is done by replacing every occourence of the type double in the function with the type #RealForward.
The overloading implementation looks now like
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

      stc::cout << "f(4.0) = " << y << std::endl;

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
For the function \f$y = f(x)\f$  with the input variables \f$x \in \mathbb{R}^n\f$ and the output variables \f$y \in \mathbb{R}^m\f$.
The forward mode of AD computes the equation
\f[
  \dot y = \frac{\partial f}{\partial x}(x) \cdot \dot x \quad .
\f]
In our example the dimensions of the vectorspaces \f$\mathbb{R}^n\f$ and \f$\mathbb{R}^n\f$ is one.
The formulation contains the partial derivative of \f$f\f$ with respect to the input variable \f$x\f$.
At the start of the tutorial this partial derivative for this simple function was derived.
We see now, that AD will multiply this partial derivative with the value of \f$\dot x\f$ and the result
of the operation will be in \f$\dot y\f$.
In the implementation for the overloaded function we have only the definition of x as the type RealForward.
But this type contains the value for the primal and the derivative which is in this case \f$\dot x\f$.
The primal value can be accessed with the functions getValue(), setValue() and value().
For the gradient the same functions getGradient(), setGradient() and gradient() exist.
The value of \f$\dot x\f$ can be set with the call:
~~~~{.cpp}
      x.setGradient(1.0);
~~~~
The value of \f$\dot y\f$ can be retrived with the call:
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

The result for the primal computation will be 64 and for the derivative it will be 48.

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
