/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2019 Chair for Scientific Computing (SciComp), TU Kaiserslautern
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

#include "../output.hpp"

void evalTest(std::vector<NUMBER>& x, std::vector<NUMBER>& y) {
  func(x.data(), y.data());
}


int main(int nargs, char** args) {
  (void)nargs;
  (void)args;

  int evalPoints = getEvalPointsCount();
  int inputs = getInputCount();
  int outputs = getOutputCount();

  std::vector<double> x(inputs);

  auto handle = codi::EvaluationHelper::template create<NUMBER>(&evalTest, outputs, inputs);
  auto jac = codi::EvaluationHelper::createJacobian(outputs, inputs);

  for(int curPoint = 0; curPoint < evalPoints; ++curPoint) {
    std::cout << "Point " << curPoint << " : {";

    for(int i = 0; i < inputs; ++i) {
      if(i != 0) {
        std::cout << ", ";
      }
      x[i] = getEvalPoint(curPoint, i);
      std::cout << x[i];
    }
    std::cout << "}\n";

    codi::EvaluationHelper::evalHandleJacobian(handle, x, jac);

    // evaluate a second time to force at least one tape reset.
    codi::EvaluationHelper::evalHandleJacobian(handle, x, jac);

    writeOutputJacobian(jac);
  }
}
