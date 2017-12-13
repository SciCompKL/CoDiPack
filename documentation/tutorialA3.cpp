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
