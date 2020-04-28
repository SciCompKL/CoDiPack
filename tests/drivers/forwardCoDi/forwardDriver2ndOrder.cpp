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

  int evalPoints = getEvalPointsCount();
  int inputs = getInputCount();
  int outputs = getOutputCount();
  NUMBER* x = new NUMBER[inputs];
  NUMBER* y = new NUMBER[outputs];

  codi::Hessian<std::vector<double>> hes(outputs, inputs);

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

    for(int curIn1st = 0; curIn1st < inputs; ++curIn1st) {
      x[curIn1st].gradient().value() = 1.0;

      for(int curIn2nd = 0; curIn2nd < inputs; ++curIn2nd) {
        x[curIn2nd].value().gradient() = 1.0;

        func(x, y);

        for(int curOut = 0; curOut < outputs; ++curOut) {
          hes(curOut, curIn1st, curIn2nd) =  y[curOut].getGradient().getGradient();
        }

        x[curIn2nd].value().setGradient(Real::Real());
      }

      x[curIn1st].setGradient(Gradient());
    }

    writeOutputHessian(hes);
  }

  delete [] x;
  delete [] y;
}
