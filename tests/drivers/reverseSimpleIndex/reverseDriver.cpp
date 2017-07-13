/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2017 Chair for Scientific Computing (SciComp), TU Kaiserslautern
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
 * Authors: Max Sagebaum, Tim Albring, (SciComp, TU Kaiserslautern)
 */

#include <toolDefines.h>

#include <iostream>
#include <vector>

int main(int nargs, char** args) {
  int evalPoints = getEvalPointsCount();
  int inputs = getInputCount();
  int outputs = getOutputCount();
  NUMBER* x = new NUMBER[inputs];
  NUMBER* y = new NUMBER[outputs];

  NUMBER::TapeType& tape = NUMBER::getGlobalTape();
  tape.resize(1000, 1000);
  tape.setExternalFunctionChunkSize(1000);

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

    std::vector<std::vector<double> > jac(outputs);
    for(int curOut = 0; curOut < outputs; ++curOut) {
      tape.setActive();
      for(int i = 0; i < inputs; ++i) {
        tape.registerInput(x[i]);
      }

      func(x, y);

      for(int i = 0; i < outputs; ++i) {
        tape.registerOutput(y[i]);
      }

      tape.setPassive();

      for(int i = 0; i < outputs; ++i) {
        y[i].setGradient(i == curOut ? 1.0:0.0);
      }

      tape.evaluate();

      for(int curIn = 0; curIn < inputs; ++curIn) {
        jac[curOut].push_back(x[curIn].getGradient());
      }

      tape.reset();
    }

    for(int curIn = 0; curIn < inputs; ++curIn) {
      for(int curOut = 0; curOut < outputs; ++curOut) {
        std::cout << curIn << " " << curOut << " " << jac[curOut][curIn] << std::endl;
      }
    }
  }
}
