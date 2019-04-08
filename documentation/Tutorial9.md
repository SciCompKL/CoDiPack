Tutorial 9: Gradient evaluations for different program configurations {#Tutorial9}
============

This tutorial describes how different parts of a program can be differentiated one after another or
how the same part in a program with different configurations can be handled.
The main difficulty with CoDiPack is to ensure that leftovers from previous recordings are properly cleaned before recording again.

The differentiated function for this tutorial is quite simple, since we want to focus on the CoDiPack specific aspects:
~~~~{.cpp}
    double func(const double& a, const double& b) {
      return a * b;
    }
~~~~

As the function is quite simple the Jacobian of the function can be computed by hand and is
\f[
  \frac{\d f}{\d a}(a,b) = b \eqdot
\f]
\f[
  \frac{\d f}{\d b}(a,b) = a \eqdot
\f]
In order to demonstrate the cleaning up of the previous runs, we want to differentiate the function once with respect to \f$a\f$
and once with respect to \f$b\f$.

The approach from the previous tutorials would yield the code:
The recording and evaluation procedure from Tutorial 2 is:
~~~~{.cpp}
      codi::RealReverse::TapeType& tape = codi::RealReverse::getGlobalTape();
      codi::RealReverse a = 4.0;
      codi::RealReverse b = 3.0;

      // record with respect to a
      tape.setActive();

      tape.registerInput(a);
      codi::RealReverse y = func(a, b);
      tape.registerOutput(y);

      tape.setPassive();
      y.setGradient(1.0);
      tape.evaluate();

      std::cout << "f(4.0, 3.0) = " << y << std::endl;
      std::cout << "df/da(4.0, 3.0) = " << a.getGradient() << std::endl;

      // record with respect to b
      tape.setActive();

      tape.registerInput(b);
      codi::RealReverse y = func(a, b);
      tape.registerOutput(y);

      tape.setPassive();
      y.setGradient(1.0);
      tape.evaluate();

      std::cout << "f(4.0, 3.0) = " << y << std::endl;
      std::cout << "df/da(4.0, 3.0) = " << a.getGradient() << std::endl;
~~~~

The output from this procedure is:
~~~~
f(4.0, 3.0) = 12
df/da(4.0, 3.0) = 3
df/db(4.0, 3.0) = 0
f(4.0, 3.0) = 12
df/da(4.0, 3.0) = 7
df/db(4.0, 3.0) = 7
~~~~
The first derivative with respect to a is correct but the second derivative with respect to b is not correct.
Also the gradient from \f$a\f$ has still a value which should not be the case.

The problem stems from the internal representation of AD in the implementation of CoDiPack.
Each variable is associated with an index and this index is used to look up the gradient.
The zero index has special meaning for CoDiPack. Each variable with a zero index does not
depend on any input variables which is also the default. The 'registerInput' function will for example
assign a new index to the given variable.

The wrong gradient from the above example comes from the leftovers from the first recording.
\f$a\f$ has still an associated index from the first gradient recording. The 'reset' call on the tape
will also reset the index handling of the tape and therefore \f$b\f$ will receive the same index as \$fa\$f since there are
both the first variable with are defined as inputs. The second call of 'func' sees now \f$a\f$ and \f$b\f$ as active variables.
Since both have the same index the reverse evaluation in CoDiPack treats them as the same adjoint variable.

The solution to the problem is to clear the AD specific information from \f$a\f$. The function for this is 'deactivateValue' from
the reverse tape interface.

The changed code will now look like:
~~~~{.cpp}
  codi::RealReverse::TapeType& tape = codi::RealReverse::getGlobalTape();
  codi::RealReverse a = 4.0;
  codi::RealReverse b = 3.0;

  // record with respect to a
  tape.setActive();

  tape.registerInput(a);
  codi::RealReverse y = func(a, b);
  tape.registerOutput(y);

  tape.setPassive();
  y.setGradient(1.0);
  tape.evaluate();

  std::cout << "f(4.0, 3.0) = " << y << std::endl;
  std::cout << "df/da(4.0, 3.0) = " << a.getGradient() << std::endl;
  std::cout << "df/db(4.0, 3.0) = " << b.getGradient() << std::endl;

  tape.deactivateValue(a);

  // record with respect to b
  tape.reset();
  tape.setActive();

  tape.registerInput(b);
  y = func(a, b);
  tape.registerOutput(y);

  tape.setPassive();
  y.setGradient(1.0);
  tape.evaluate();

  std::cout << "f(4.0, 3.0) = " << y << std::endl;
  std::cout << "df/da(4.0, 3.0) = " << a.getGradient() << std::endl;
  std::cout << "df/db(4.0, 3.0) = " << b.getGradient() << std::endl;
~~~~
and the result is:
~~~~
f(4.0, 3.0) = 12
df/da(4.0, 3.0) = 3
df/db(4.0, 3.0) = 0
f(4.0, 3.0) = 12
df/da(4.0, 3.0) = 0
df/db(4.0, 3.0) = 4
~~~~

The general in applications which are differentiated with CoDiPack is, that all values from previous gradient computations should either be
registered as inputs, overwritten in the computation or deactivated with the above method. The above example shows the behavior of the default
codi::RealReverse type. All types without an 'Index' in there name will behave the same. That is, if values from a previous gradient computation
are used then the derivative results will be wrong.

The behavior for the codi::RealReverseIndex type and all types with an 'Index' in there names is different. The correctness of the gradient for
these types is not affected by old values from previous gradient computations. Nevertheless, the old values are still seen as value that
depend on the inputs and are therefore treated by CoDiPack. For large application this can impact the performance of the gradient computation.

The full code of the example is:
~~~~{.cpp}
#include <codi.hpp>
#include <iostream>

codi::RealReverse func(const codi::RealReverse& a, const codi::RealReverse& b) {
  return a * b;
}

void call(bool clear, bool stats) {
  codi::RealReverse::TapeType& tape = codi::RealReverse::getGlobalTape();
  codi::RealReverse a = 4.0;
  codi::RealReverse b = 3.0;

  // record with respect to a
  tape.setActive();

  tape.registerInput(a);
  codi::RealReverse y = func(a, b);
  tape.registerOutput(y);

  tape.setPassive();
  y.setGradient(1.0);
  tape.evaluate();

  std::cout << "f(4.0, 3.0) = " << y << std::endl;
  std::cout << "df/da(4.0, 3.0) = " << a.getGradient() << std::endl;
  std::cout << "df/db(4.0, 3.0) = " << b.getGradient() << std::endl;

  if(stats) {
    tape.printStatistics();
  }
  if(clear) {
    tape.deactivateValue(a);
  }

  // record with respect to b
  tape.reset();
  tape.setActive();

  tape.registerInput(b);
  y = func(a, b);
  tape.registerOutput(y);

  tape.setPassive();
  y.setGradient(1.0);
  tape.evaluate();

  std::cout << "f(4.0, 3.0) = " << y << std::endl;
  std::cout << "df/da(4.0, 3.0) = " << a.getGradient() << std::endl;
  std::cout << "df/db(4.0, 3.0) = " << b.getGradient() << std::endl;

  if(stats) {
    tape.printStatistics();
  }
}
int main(int nargs, char** args) {

  std::cout << "Recording tapes without clear:" << std::endl;
  call(false, false);

  std::cout << std::endl;
  std::cout << "Recording tapes with clear:" << std::endl;
  call(true, false);

  return 0;
}
~~~~
