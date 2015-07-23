Tutorial 2: Reverse mode      {#Tutorial2}
============


The function we want to differentiate is implemented as
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
  \tag{T2.1}
  \frac{\d f}{\d x}(x) = 3 * x * x = 3 * x^2 \eqdot
\f]
The derivative computation with the reverse mode of CoDiPack needs now seven steps:
  - An implementation of the function with the CoDiPack reverse type [RealReverse](@ref codi::RealReverse).
  - A definition of the input variables.
  - Record the function evaluation on the tape.
  - A definition of the output variables.
  - The direction of the derivative has to be set.
  - The tape recorded by CoDiPack has to be evaluated.
  - The result of the derivative has to be received.

The first step is the implementation of the function with the reverse type of CoDiPack.
This is done by replacing every occurrence of the type double in the function with the type RealReverse.
The overloading implementation is
~~~~{.cpp}
    codi::RealReverse func(const codi::RealReverse& x) {
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

The overloaded function can only be called if we exchange the double types with the RealReverse type.
~~~~{.cpp}
    int main(int nargs, char** args) {
      codi::RealReverse x = 4.0;

      codi::RealReverse y = func(x);

      std::cout << "f(4.0) = " << y << std::endl;

      return 0;
    }
~~~~

Now the implementation of the function with the RealReverse type is completed.
The next step will register the input value for the function on the tape structure for the RealReverse type.
The tape can be accessed with the getGlobalTape function.
The common use case will store the tape in a reference and then perform operations on the tape.
The code which gets a reference to the tape, activates it and registers the input looks as:
~~~~{.cpp}
   codi::RealReverse::TapeType& tape = codi::RealReverse::getGlobalTape();
   tape.setActive();
   tape.registerInput(x);

   codi::RealReverse y = func(x);
   tape.registerOutput(y);
   tape.setPassive();
~~~~
The last three lines in the code block call the function evaluation which records the tape for \f$f\f$ and
register the output variable.
Afterwards the tape is deactivated and no further statements will be recorded.
The next step requires the definition of the gradient direction during the tape evaluation.
In order to get a better feeling what this means, the reverse AD formulation is repeated here.
For the function \f$y = f(x)\f$  with the input variables \f$x \in \R^n\f$ and the output variables \f$y \in \R^m\f$.
The reverse mode of AD computes the equation
\f[
  \bar x = \frac{\d f}{\d x}^T(x) \cdot \bar y \quad .
\f]
In our example the dimensions of the vector spaces \f$\R^n\f$ and \f$\R^m\f$ is one.
The formulation contains the partial derivative of \f$f\f$ with respect to the input variable \f$x\f$,
which was derived in equation (T2.1).
AD will multiply this partial derivative with the value of \f$\bar y\f$ and the result of the operation will be stored in \f$\bar x\f$.
By setting \f$\bar y = 1\f$ the result will be the gradient of \f$f\f$ at the position \f$x\f$.
One important aspect is the reversal of the input and output values.
For the original function \f$x\f$ is the input value and \f$y\f$ is the output value.
In the reverse mode the roles are switched, \f$\bar y\f$ is an input value and \f$\bar x\f$ is an output value.

The RealReverse type provides functions to access the primal value and the derivative value which is denoted by the bar above the name of the value.
The primal value can be accessed with the functions getValue(), setValue() and value().
For the gradient the same functions getGradient(), setGradient() and gradient() exist.
The value of \f$\bar y\f$ can be set with the call:
~~~~{.cpp}
      y.setGradient(1.0);
~~~~
The value of \f$\bar x\f$ can be retrieved after the tape has been evaluated.
The calls for this are
~~~~{.cpp}
      tape.evaluate();
      double x_bar = x.getGradient();
~~~~

The implementation of the main method is now extended by these calls:
~~~~{.cpp}
    int main(int nargs, char** args) {
      codi::RealReverse x = 4.0;

      codi::RealReverse::TapeType& tape = codi::RealReverse::getGlobalTape();
      tape.setActive();

      tape.registerInput(x);
      codi::RealReverse y = func(x);
      tape.registerOutput(y);

      tape.setPassive();
      y.setGradient(1.0);
      tape.evaluate();

      std::cout << "f(4.0) = " << y << std::endl;
      std::cout << "df/dx(4.0) = " << x.getGradient() << std::endl;

      return 0;
    }
~~~~

The result for the primal computation will be 64 and the derivative will be 48.

The full code of the example is now:
~~~~{.cpp}
    #include <codi.hpp>
    #include <iostream>

    codi::RealReverse func(const codi::RealReverse& x) {
      return x * x * x;
    }

    int main(int nargs, char** args) {
      codi::RealReverse x = 4.0;

      codi::RealReverse::TapeType& tape = codi::RealReverse::getGlobalTape();
      tape.setActive();

      tape.registerInput(x);
      codi::RealReverse y = func(x);
      tape.registerOutput(y);

      tape.setPassive();
      y.setGradient(1.0);
      tape.evaluate();

      std::cout << "f(4.0) = " << y << std::endl;
      std::cout << "df/dx(4.0) = " << x.getGradient() << std::endl;

      return 0;
    }

~~~~
