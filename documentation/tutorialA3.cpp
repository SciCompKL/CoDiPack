/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2021 Chair for Scientific Computing (SciComp), TU Kaiserslautern
 * Homepage: http://www.scicomp.uni-kl.de
 * Contact:  Prof. Nicolas R. Gauger (codi@scicomp.uni-kl.de)
 *
 * Lead developers: Max Sagebaum, Tim Albring (SciComp, TU Kaiserslautern)
 *
 * This file is part of CoDiPack (http://www.scicomp.uni-kl.de/software/codi).
 *
 * CoDiPack is free software: you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * CoDiPack is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty
 * of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 * See the GNU General Public License for more details.
 * You should have received a copy of the GNU
 * General Public License along with CoDiPack.
 * If not, see <http://www.gnu.org/licenses/>.
 *
 * Authors:
 *  - SciComp, TU Kaiserslautern:
 *     Max Sagebaum
 *     Tim Albring
 *     Johannes Blühdorn
 */
#include <codi.hpp>

#include <iostream>
#include <algorithm>

// Evaluates w = (1 y y^2 ... y^(n-1)) A (1 x x^2 ... x^(n-1))^T
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
