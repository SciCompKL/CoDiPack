Tutorial A1: External functions {#TutorialA1}
============

If AD is applied on a large code most often this code base uses other
libraries for the computation. The most common example is the use of a
linear algebra library like [Eigen](eigen.tuxfamily.org/index.php) or [LAPACK](http://www.netlib.org/lapack/).

Lets say we want to solve the system
\f[
  x = A^{-1}b \eqdot
\f]
for \f$x\f$. The forward AD mode for this operation is
\f[
  \dot x = A^{-1}(\dot b - \dot Ax)
\f]
and the reverse AD mode for this operation is
\f[
  \tag{TA1.1}
  \begin{aligned}
      s = & A^{-T}\bar x\\
      \bar W \aeq & -s \cdot x^T \\
      \bar b \aeq & s \\
      \bar x = & 0 \eqdot
  \end{aligned}
\f]
Usually the call of a library routine to LAPACK is encapsulated in a function
call. For our example we assume that the function is defined as
~~~~{.cpp}
    template<typename Real>
    void solve2(const Real* A, const Real* b, const Real* x) {

      // A = a[0] a[1]  A^-1 = 1/det *  a[3] -a[1]
      //     a[2] a[3]                 -a[2]  a[0]
      Real det = A[0] * A[3] - A[1] * A[2];

      x[0] = (A[3] * b[0] - A[1] * b[1]) / det;
      x[1] = (-A[2] * b[0] + A[0] * b[1]) / det;
    }
~~~~
for the solution of a system with the dimension 2. We do not use any
library here to keep the example self contained.

A simple code that uses this function could look like:
~~~~{.cpp}

    double u = 3.0;

    double A[4] = {u * 1.0, 0.5,  0.25, u * -1.0};
    double b[2] = {u * 10.0, u * 20.0};

    double x[2];

    solve2(A, b, x);

    double w = sqrt(x[0] * x[0] + x[1] * x[1]);

    std::cout << "Solution w: " << w << std::endl;
~~~~

The derivative for this simple example can be computed with CoDiPack by recording the tape:
~~~~{.cpp}
    codi::RealReverse u = 3.0;

    codi::RealReverse::TapeType& tape = codi::RealReverse::getGlobalTape();
    tape.setActive();
    tape.registerInput(u);

    codi::RealReverse A[4] = { u * 1.0, 0.5,  0.25, u * -1.0};
    codi::RealReverse b[2] = {u * 10.0, u * 20.0};

    codi::RealReverse x[2];

    solve2(A, b, x);

    codi::RealReverse w = sqrt(x[0] * x[0] + x[1] * x[1]);

    tape.registerOutput(w);

    tape.setPassive();
    w.setGradient(1);

    tape.evaluate();

    std::cout << "Solution w: " << w << std::endl;
    std::cout << "Adjoint u: " << u.getGradient() << std::endl;
~~~~

Now we assume that `solve2` is implemented with an external library and
we can not differentiate it with CoDiPack. Therefore we need to use the
external function helper (codi::ExternalFunctionHelper) from CoDiPack. It
provides a standardized interface for defining which are the inputs to
an external function, the outputs and to provide additional data. The user
then needs to implement evaluation functions such that the reverse mode
can be evaluated.

The external function helper assumes that the function computes
\f[
    y = f(x) \eqdot
\f]
The reverse AD mode for this function is
\f[
    \begin{aligned}
        \bar x \aeq & \frac{\d f}{\d x}(x) \bar y \\
        \bar y = & 0
    \end{aligned}
\f]
The user needs to provide functions that compute \f$f\f$ and the reverse
mode of \f$f\f$.

The most common use case for the codi::ExternalFunctionHelper involves
the following steps:
 - Register the inputs
 - Register the outputs
 - Add the user data
 - Call the primal evaluation
 - Register the external function helper

 The first two steps are quite simple and since we do not need any user
 data we skip the third step.
~~~~{.cpp}
    codi::ExternalFunctionHelper<codi::RealReverse> eh;

    for(int i = 0; i < 4; ++i) {
      eh.addInput(A[i]);
    }
    for(int i = 0; i < 2; ++i) {
      eh.addInput(b[i]);
    }

    for(int i = 0; i < 2; ++i) {
      eh.addOutput(x[i]);
    }
~~~~

 Now the helper know about the input and output values for the code section
 that requires the special treatment. In the next step a function needs
 to be specified that can evaluate the primal computation without the
 CoDiPack types.
~~~~{.cpp}
  void solve2_primal(const codi::RealReverse::Real* x, size_t m, codi::RealReverse::Real* y, size_t n, codi::DataStore* d) {
    CODI_UNUSED(m);
    CODI_UNUSED(n);
    CODI_UNUSED(d);

    solve2(&x[0], &x[4], y);
  }
~~~~
The function `solve2_primal` will be called by the external function
helper. The helper extracts the primal values from all the input values
into one large array and provides it to the call.
This implementation calls just `solve2`. Since
all inputs are stored in one big array, we need to access the correct
array pointers in the call. Also the variable naming can be a little
bit confusing since the external function helper assumes the naming
\f$ y = f(x) \f$ and the linear system solve is defined as \f$ x = A^{-1} \f$.
In the first case \f$x\f$ describes the input values and in the second case
\f$x\f$ describes the output values.

The primal function is then called via the external function helper:
~~~~{.cpp}
  eh.callPrimalFunc(solve2_primal);
~~~~

The final step adds the external function helper to the tape:
~~~~{.cpp}
  eh.addToTape(solve2_rev);
~~~~
The function `solve2_rev` contains the logic for the reverse evaluation and
still needs to be implemented. The algorithm for this was already shown in
equation (TA1.1). The implementation is:
~~~~{.cpp}
void solve2_rev(const codi::RealReverse::Real* x, codi::RealReverse::Real* x_b, size_t m, const codi::RealReverse::Real* y, const codi::RealReverse::Real* y_b, size_t n, codi::DataStore* d) {
  CODI_UNUSED(m);
  CODI_UNUSED(n);
  CODI_UNUSED(y);
  CODI_UNUSED(d);

  codi::RealReverse::Real ATrans[4] = {x[0], x[2], x[1], x[3]};

  codi::RealReverse::Real s[2];
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
~~~~
It is important to note that we did not perform the updates (`+=`) on the
adjoint input variables nor the reset of the adjoint output variables.
The external function helper takes care of these things, so they are not
necessary.

The full code with the external function helper is:
~~~~{.cpp}
  codi::RealReverse u = 3.0;

  codi::RealReverse::TapeType& tape = codi::RealReverse::getGlobalTape();
  tape.setActive();
  tape.registerInput(u);

  codi::RealReverse A[4] = { u * 1.0, 0.5,  0.25, u * -1.0};
  codi::RealReverse b[2] = {u * 10.0, u * 20.0};

  codi::RealReverse x[2];

  // external function helper start
  codi::ExternalFunctionHelper<codi::RealReverse> eh;
  for(int i = 0; i < 4; ++i) {
    eh.addInput(A[i]);
  }
  for(int i = 0; i < 2; ++i) {
    eh.addInput(b[i]);
  }

  for(int i = 0; i < 2; ++i) {
    eh.addOutput(x[i]);
  }

  eh.callPrimalFunc(solve2_primal);
  eh.addToTape(solve2_rev);
  // external function helper end

  codi::RealReverse w = sqrt(x[0] * x[0] + x[1] * x[1]);

  tape.registerOutput(w);

  tape.setPassive();
  w.setGradient(1);

  tape.evaluate();

  std::cout << "Solution w: " << w << std::endl;
  std::cout << "Adjoint u: " << u.getGradient() << std::endl;
~~~~

In this code the `solve2` function is never called directly. The primal
solution is computed through a call of `solve2_primal` and the adjoint
solution is computed by a call of `solve2_rev` during the reverse
evaluation.

This is the first use case for the external function helper. The second
use case is less common and applied if the code section, in this case `solve2`,
can be evaluated with the CoDiPack types but a special adjoint implementation
is available. The basic concept of this operation mode is that the tape
is set to passive during the evaluation of `solve2` and the user provides
the function for the adjoint computation.
The steps for this mode are:
 - Register the inputs
 - Call the passive evaluation
 - Register the outputs
 - Add the user data
 - Register the external function helper

Since the changes in the implementation are only minor, the full code for
the external function helper is provided directly:
~~~~{.cpp}
   codi::ExternalFunctionHelper<codi::RealReverse> eh;
   for(int i = 0; i < 4; ++i) {
     eh.addInput(A[i]);
   }
   for(int i = 0; i < 2; ++i) {
     eh.addInput(b[i]);
   }

   eh.callPassiveFunc(solve2<codi::RealReverse>, A, b, x);

   for(int i = 0; i < 2; ++i) {
     eh.addOutput(x[i]);
   }

   eh.addToTape(solve2_rev);
~~~~
The major differences are that the outputs are registered after the
passive function call and there is no need to provide the primal
evaluation function `solve2_primal`.

The full code for the tutorial is:
~~~~{.cpp}
#include <codi.hpp>

#include <iostream>

template<typename Real>
void solve2(const Real* A, const Real* b, Real* x) {

  // A = a[0] a[1]  A^-1 = 1/det *  a[3] -a[1]
  //     a[2] a[3]                 -a[2]  a[0]
  Real det = A[0] * A[3] - A[1] * A[2];

  x[0] = (A[3] * b[0] - A[1] * b[1]) / det;
  x[1] = (-A[2] * b[0] + A[0] * b[1]) / det;
}

void solve2_primal(const codi::RealReverse::Real* x, size_t m, codi::RealReverse::Real* y, size_t n, codi::DataStore* d) {
  CODI_UNUSED(m);
  CODI_UNUSED(n);
  CODI_UNUSED(d);

  solve2(&x[0], &x[4], y);
}

void solve2_rev(const codi::RealReverse::Real* x, codi::RealReverse::Real* x_b, size_t m, const codi::RealReverse::Real* y, const codi::RealReverse::Real* y_b, size_t n, codi::DataStore* d) {
  CODI_UNUSED(m);
  CODI_UNUSED(n);
  CODI_UNUSED(y);
  CODI_UNUSED(d);

  codi::RealReverse::Real ATrans[4] = {x[0], x[2], x[1], x[3]};

  codi::RealReverse::Real s[2];
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

void primal() {
  std::cout << "double:" << std::endl;

  double u = 3.0;

  double A[4] = {u * 1.0, 0.5,  0.25, u * -1.0};
  double b[2] = {u * 10.0, u * 20.0};

  double x[2];

  solve2(A, b, x);

  double w = sqrt(x[0] * x[0] + x[1] * x[1]);

  std::cout << "Solution w: " << w << std::endl;
}

void derivative() {
  std::cout << "codi::RealReverse:" << std::endl;

  codi::RealReverse u = 3.0;

  codi::RealReverse::TapeType& tape = codi::RealReverse::getGlobalTape();
  tape.setActive();
  tape.registerInput(u);

  codi::RealReverse A[4] = { u * 1.0, 0.5,  0.25, u * -1.0};
  codi::RealReverse b[2] = {u * 10.0, u * 20.0};

  codi::RealReverse x[2];

  solve2(A, b, x);

  codi::RealReverse w = sqrt(x[0] * x[0] + x[1] * x[1]);

  tape.registerOutput(w);

  tape.setPassive();
  w.setGradient(1);

  tape.evaluate();

  std::cout << "Solution w: " << w << std::endl;
  std::cout << "Adjoint u: " << u.getGradient() << std::endl;
}

void externalFunction() {
  std::cout << "codi::RealReverse(External function):" << std::endl;

  codi::RealReverse u = 3.0;

  codi::RealReverse::TapeType& tape = codi::RealReverse::getGlobalTape();
  tape.setActive();
  tape.registerInput(u);

  codi::RealReverse A[4] = { u * 1.0, 0.5,  0.25, u * -1.0};
  codi::RealReverse b[2] = {u * 10.0, u * 20.0};

  codi::RealReverse x[2];

  // external function helper start
  codi::ExternalFunctionHelper<codi::RealReverse> eh;
  for(int i = 0; i < 4; ++i) {
    eh.addInput(A[i]);
  }
  for(int i = 0; i < 2; ++i) {
    eh.addInput(b[i]);
  }

  for(int i = 0; i < 2; ++i) {
    eh.addOutput(x[i]);
  }

  eh.callPrimalFunc(solve2_primal);
  eh.addToTape(solve2_rev);
  // external function helper end

  codi::RealReverse w = sqrt(x[0] * x[0] + x[1] * x[1]);

  tape.registerOutput(w);

  tape.setPassive();
  w.setGradient(1);

  tape.evaluate();

  std::cout << "Solution w: " << w << std::endl;
  std::cout << "Adjoint u: " << u.getGradient() << std::endl;
}

void externalFunctionPassive() {
  std::cout << "codi::RealReverse(External function passive):" << std::endl;

  codi::RealReverse u = 3.0;

  codi::RealReverse::TapeType& tape = codi::RealReverse::getGlobalTape();
  tape.setActive();
  tape.registerInput(u);

  codi::RealReverse A[4] = { u * 1.0, 0.5,  0.25, u * -1.0};
  codi::RealReverse b[2] = {u * 10.0, u * 20.0};

  codi::RealReverse x[2];

  // external function helper start
  codi::ExternalFunctionHelper<codi::RealReverse> eh;
  for(int i = 0; i < 4; ++i) {
    eh.addInput(A[i]);
  }
  for(int i = 0; i < 2; ++i) {
    eh.addInput(b[i]);
  }

  eh.callPassiveFunc(solve2<codi::RealReverse>, A, b, x);

  for(int i = 0; i < 2; ++i) {
    eh.addOutput(x[i]);
  }

  eh.addToTape(solve2_rev);
  // external function helper end

  codi::RealReverse w = sqrt(x[0] * x[0] + x[1] * x[1]);

  tape.registerOutput(w);

  tape.setPassive();
  w.setGradient(1);

  tape.evaluate();

  std::cout << "Solution w: " << w << std::endl;
  std::cout << "Adjoint u: " << u.getGradient() << std::endl;
}

int main(int nargs, char** args) {

  primal();
  derivative();
  externalFunction();
  externalFunctionPassive();

  return 0;
}
~~~~
