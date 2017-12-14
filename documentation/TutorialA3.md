Tutorial A3: Manual statement storing {#TutorialA3}
============

CoDiPack stores for every statement in the in the code some data such that
it can perform the reverse evaluation. Sometimes a simple function evaluation
consists of multiple statements but has only few input values and one
output value.

For the function
\f[
    \tag{TA3.1}
    y = f(x)
\f]
with \f$x \in \R^n\f$ and \f$y \in \R\f$, the reverse AD mode is
\f[
    \tag{TA3.2}
    \begin{aligned}
        \bar x \aeq & \frac{\d f}{\d x}(x)^T \bar y \\
        \bar y = & 0
    \end{aligned}
\f]
A simple function could be the evaluation of a polynomial:
~~~~{.cpp}
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
~~~~

The derivative with respect to `x` and `y` can be programmed in a very simple
way and we could compute the Jacobian \f$\frac{\d f}{\d (x,y)}\f$. This
Jacobian can then be stored on the tape instead of the operations that
CoDiPack would usually store. The codi::StatementPushHelper can be used
to perform this task.

The context for the evaluation of the polynomial is:
~~~~{.cpp}
    codi::RealReverse u = 3.0;

    codi::RealReverse::TapeType& tape = codi::RealReverse::getGlobalTape();
    tape.setActive();
    tape.registerInput(u);

    double A[9] = { 1.0,    0.5,   0.25,
                    0.0,    1.0,   0.75,
                    0.25,   0.0,   1.0};
    codi::RealReverse x = cos(u);
    codi::RealReverse y = sin(u);

    codi::RealReverse o = poly2D(x, y, A, 3);


    codi::RealReverse w = exp(o * o);

    tape.registerOutput(w);

    tape.setPassive();
    w.setGradient(1);

    tape.evaluate();

    tape.printStatistics();

    std::cout << "Solution w: " << w << std::endl;
    std::cout << "Adjoint u: " << u.getGradient() << std::endl;
~~~~
The complete recording of this section uses 704 bytes.

In order to compute the the derivatives of the Jacobi, the derivative computation
needs to be available:
~~~~{.cpp}
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
~~~~

The steps for storing the Jacobian in (TA3.2) for the equation (TA3.1)
with the codi::StatementPushHelper are the following:
 - Compute y in a passive context
 - Compute \f$\frac{\d f}{\d x}\f$ in a passive context
 - Add the Jacobi to the tape.

The passive context is important at this point. Otherwise CoDiPack would
still record the operations if the CoDiPack type is used.

The first example will use the double type for the evaluation.
The computation of the polynomial value and the Jacobian is quite simple:
~~~~{.cpp}
  double o_p = poly2D(x.getValue(), y.getValue(), A, 3);
  double jac[2];
  jac[0] = poly2D_dx(x.getValue(), y.getValue(), A, 3);
  jac[1] = poly2D_dy(x.getValue(), y.getValue(), A, 3);
~~~~

The next step is now the storing of the Jacobian on the tape:
~~~~{.cpp}
  codi::StatementPushHelper<codi::RealReverse> ph;

  ph.startPushStatement();
  ph.pushArgument(x, jac[0]);
  ph.pushArgument(y, jac[1]);
  ph.endPushStatement(o, o_p);
~~~~

The number of arguments for one expression is capped at 255, otherwise the
statement push helper will ensure that the variables are handled in the
correct way.

The full code for the memory optimized example is:
~~~~{.cpp}
    codi::RealReverse u = 3.0;

    codi::RealReverse::TapeType& tape = codi::RealReverse::getGlobalTape();
    tape.setActive();
    tape.registerInput(u);

    double A[9] = { 1.0,    0.5,   0.25,
                    0.0,    1.0,   0.75,
                    0.25,   0.0,   1.0};
    codi::RealReverse x = cos(u);
    codi::RealReverse y = sin(u);

    // start manual statement push
    codi::RealReverse o;
    double o_p = poly2D(x.getValue(), y.getValue(), A, 3);
    double jac[2];
    jac[0] = poly2D_dx(x.getValue(), y.getValue(), A, 3);
    jac[1] = poly2D_dy(x.getValue(), y.getValue(), A, 3);

    codi::StatementPushHelper<codi::RealReverse> ph;

    ph.startPushStatement();
    ph.pushArgument(x, jac[0]);
    ph.pushArgument(y, jac[1]);
    ph.endPushStatement(o, o_p);
    //end manual statement push

    codi::RealReverse w = exp(o * o);

    tape.registerOutput(w);

    tape.setPassive();
    w.setGradient(1);

    tape.evaluate();

    tape.printStatistics();

    std::cout << "Solution w: " << w << std::endl;
    std::cout << "Adjoint u: " << u.getGradient() << std::endl;
~~~~
The memory used by this code is only 146 bytes. This includes also the
computations that are evaluated before and after the polynomial evaluation.
If only the memory for the polynomial evaluation is counted then the
original evaluation used 591 bytes and the optimized evaluation uses 33
bytes, which is a factor of 18.

The same result could also be archived with the external functions from
[tutorial A1](@ref TutorialA1) but there additional memory is required for
the external function. This would diminish the improvement of the memory
optimization.

Sometimes it is not possible to evaluate the primal function and/or the
derivative computation with double values. If they need to be evaluated
with CoDiPack type the tape can be set to passive during the computation
of the primal and gradient values. Afterwards, the tape needs to be set to active
again. In order for the statement push helper to work the tape needs to be
active during startPushStatement and endPushStatement methods.
For the current example the code would look like:
~~~~{.cpp}
    codi::RealReverse::getGlobalTape().setPassive();
    codi::RealReverse o = poly2D(x, y, A, 3);
    codi::RealReverse jac[2];
    jac[0] = poly2D_dx(x, y, A, 3);
    jac[1] = poly2D_dy(x, y, A, 3);
    codi::RealReverse::getGlobalTape().setActive();

    codi::StatementPushHelper<codi::RealReverse> ph;

    ph.startPushStatement();
    ph.pushArgument(x, jac[0].getValue());
    ph.pushArgument(y, jac[1].getValue());
    ph.endPushStatement(o, o.getValue());
~~~~

There are some other convenient functions where arrays can be provided
and the statement push is evaluated in one sweep. The array example would
look like:
~~~~{.cpp}
    codi::StatementPushHelper<codi::RealReverse> ph;
    codi::RealReverse input[2] = {x, y};

    ph.pushStatement(o, o_p, input, jac, 2);
~~~~

The codi::StatementPushHelper has the following properties:
 - It can be used multiple times. Each time endPushStatement is called,
   the structure returns to a state where it can be used again.
 - The maximum number of arguments for a function is 255.

The full code for the tutorial is:
~~~~{.cpp}
#include <codi.hpp>

#include <iostream>
#include <algorithm>

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

void derivative() {
  std::cout << "codi::RealReverse:" << std::endl;

  codi::RealReverse u = 3.0;

  codi::RealReverse::TapeType& tape = codi::RealReverse::getGlobalTape();
  tape.setActive();
  tape.registerInput(u);

  double A[9] = { 1.0,    0.5,   0.25,
                  0.0,    1.0,   0.75,
                  0.25,   0.0,   1.0};
  codi::RealReverse x = cos(u);
  codi::RealReverse y = sin(u);

  codi::RealReverse o = poly2D(x, y, A, 3);


  codi::RealReverse w = exp(o * o);

  tape.registerOutput(w);

  tape.setPassive();
  w.setGradient(1);

  tape.evaluate();

  tape.printStatistics();

  std::cout << "Solution w: " << w << std::endl;
  std::cout << "Adjoint u: " << u.getGradient() << std::endl;
}

void statementPush() {
  std::cout << "codi::RealReverse(statementPush):" << std::endl;

  codi::RealReverse u = 3.0;

  codi::RealReverse::TapeType& tape = codi::RealReverse::getGlobalTape();
  tape.setActive();
  tape.registerInput(u);

  double A[9] = { 1.0,    0.5,   0.25,
                  0.0,    1.0,   0.75,
                  0.25,   0.0,   1.0};
  codi::RealReverse x = cos(u);
  codi::RealReverse y = sin(u);

  // start manual statement push
  codi::RealReverse o;
  double o_p = poly2D(x.getValue(), y.getValue(), A, 3);
  double jac[2];
  jac[0] = poly2D_dx(x.getValue(), y.getValue(), A, 3);
  jac[1] = poly2D_dy(x.getValue(), y.getValue(), A, 3);

  codi::StatementPushHelper<codi::RealReverse> ph;

  ph.startPushStatement();
  ph.pushArgument(x, jac[0]);
  ph.pushArgument(y, jac[1]);
  ph.endPushStatement(o, o_p);
  //end manual statement push

  codi::RealReverse w = exp(o * o);

  tape.registerOutput(w);

  tape.setPassive();
  w.setGradient(1);

  tape.evaluate();

  tape.printStatistics();

  std::cout << "Solution w: " << w << std::endl;
  std::cout << "Adjoint u: " << u.getGradient() << std::endl;
}

void statementPushPassive() {
  std::cout << "codi::RealReverse(statementPush Passive):" << std::endl;

  codi::RealReverse u = 3.0;

  codi::RealReverse::TapeType& tape = codi::RealReverse::getGlobalTape();
  tape.setActive();
  tape.registerInput(u);

  double A[9] = { 1.0,    0.5,   0.25,
                  0.0,    1.0,   0.75,
                  0.25,   0.0,   1.0};
  codi::RealReverse x = cos(u);
  codi::RealReverse y = sin(u);

  // start manual statement push
  codi::RealReverse::getGlobalTape().setPassive();
  codi::RealReverse o = poly2D(x, y, A, 3);
  codi::RealReverse jac[2];
  jac[0] = poly2D_dx(x, y, A, 3);
  jac[1] = poly2D_dy(x, y, A, 3);
  codi::RealReverse::getGlobalTape().setActive();

  codi::StatementPushHelper<codi::RealReverse> ph;

  ph.startPushStatement();
  ph.pushArgument(x, jac[0].getValue());
  ph.pushArgument(y, jac[1].getValue());
  ph.endPushStatement(o, o.getValue());
  //end manual statement push

  codi::RealReverse w = exp(o * o);

  tape.registerOutput(w);

  tape.setPassive();
  w.setGradient(1);

  tape.evaluate();

  tape.printStatistics();

  std::cout << "Solution w: " << w << std::endl;
  std::cout << "Adjoint u: " << u.getGradient() << std::endl;
}

void statementPushArray() {
  std::cout << "codi::RealReverse(statementPush Array):" << std::endl;

  codi::RealReverse u = 3.0;

  codi::RealReverse::TapeType& tape = codi::RealReverse::getGlobalTape();
  tape.setActive();
  tape.registerInput(u);

  double A[9] = { 1.0,    0.5,   0.25,
                  0.0,    1.0,   0.75,
                  0.25,   0.0,   1.0};
  codi::RealReverse x = cos(u);
  codi::RealReverse y = sin(u);

  // start manual statement push
  codi::RealReverse o;
  double o_p = poly2D(x.getValue(), y.getValue(), A, 3);
  double jac[2];
  jac[0] = poly2D_dx(x.getValue(), y.getValue(), A, 3);
  jac[1] = poly2D_dy(x.getValue(), y.getValue(), A, 3);

  codi::StatementPushHelper<codi::RealReverse> ph;
  codi::RealReverse input[2] = {x, y};

  ph.pushStatement(o, o_p, input, jac, 2);
  //end manual statement push

  codi::RealReverse w = exp(o * o);

  tape.registerOutput(w);

  tape.setPassive();
  w.setGradient(1);

  tape.evaluate();

  tape.printStatistics();

  std::cout << "Solution w: " << w << std::endl;
  std::cout << "Adjoint u: " << u.getGradient() << std::endl;
}

int main(int nargs, char** args) {

  derivative();
  codi::RealReverse::getGlobalTape().reset();
  statementPush();
  codi::RealReverse::getGlobalTape().reset();
  statementPushPassive();
  codi::RealReverse::getGlobalTape().reset();
  statementPushArray();

  return 0;
}
~~~~