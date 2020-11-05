/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2020 Chair for Scientific Computing (SciComp), TU Kaiserslautern
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

#include <toolDefines.h>

#include <iostream>

#include "../output.hpp"

int main(int nargs, char** args) {
  (void)nargs;
  (void)args;

  using GT = codi::GradientValueTraits<Gradient>;
  constexpr size_t gradDim = GT::getVectorSize();

  int evalPoints = getEvalPointsCount();
  int inputs = getInputCount();
  int outputs = getOutputCount();
  NUMBER* x = new NUMBER[inputs];
  NUMBER* y = new NUMBER[outputs];

  codi::Jacobian<std::vector<double>> jac(outputs, inputs);

  for(int curPoint = 0; curPoint < evalPoints; ++curPoint) {
    std::cout << "Point " << curPoint << " : {";

    for(int i = 0; i < inputs; ++i) {
      if(i != 0) {
        std::cout << ", ";
      }
      double val = getEvalPoint(curPoint, i);
      std::cout << val;

      x[i] = (NUMBER)(val);
    }
    std::cout << "}\n";

    for(int i = 0; i < outputs; ++i) {
      y[i] = 0.0;
    }

    int runs = inputs / gradDim;
    if(inputs % gradDim != 0) {
      runs += 1;
    }
    for(int curIn = 0; curIn < runs; ++curIn) {
      size_t curSize = gradDim;
      if((curIn + 1) * gradDim  > (size_t)inputs) {
        curSize = inputs % gradDim;
      }
      for(size_t curDim = 0; curDim < curSize; ++curDim) {
        GT::at(x[curIn * gradDim + curDim].gradient(), curDim) = 1.0;
      }

      for(int i = 0; i < outputs; ++i) {
        y[i].setGradient(Gradient());
      }

      func(x, y);

      for(size_t curDim = 0; curDim < curSize; ++curDim) {
        for(int curOut = 0; curOut < outputs; ++curOut) {
#if SECOND_ORDER
          jac(curOut, curIn) = y[curOut].getGradient().getValue();
#else
          jac(curOut, curIn * gradDim + curDim) = GT::at(y[curOut].getGradient(), curDim);
#endif
        }
      }

      for(size_t curDim = 0; curDim < curSize; ++curDim) {
        x[curIn * gradDim + curDim].setGradient(Gradient());
      }
    }

    writeOutputJacobian(jac);
  }

  delete [] x;
  delete [] y;
}
