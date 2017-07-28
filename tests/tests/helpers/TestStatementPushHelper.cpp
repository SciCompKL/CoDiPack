/**
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015 Chair for Scientific Computing (SciComp), TU Kaiserslautern
 * Homepage: http://www.scicomp.uni-kl.de
 * Contact:  Prof. Nicolas R. Gauger (codi@scicomp.uni-kl.de)
 *
 * Lead developers: Max Sagebaum, Tim Albring (SciComp, TU Kaiserslautern)
 *
 * This file is part of CoDiPack (http://www.scicomp.uni-kl.de/software/codi).
 *
 * CoDiPack is free software: you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation, either version 2 of the
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

#include <vector>
#include <type_traits>

IN(2)
OUT(8)
POINTS(1) = {
  {  1.0,     0.5}
};


void evalFunc(NUMBER* x, NUMBER* y) {
  y[0] = x[0];
  y[1] = x[1];
  for(int i = 0; i < 5; ++i) {
    NUMBER xTemp = y[0];
    NUMBER yTemp = y[1];

    y[0] = xTemp * xTemp - yTemp * yTemp - 0.65;
    y[1] = 2.0 * yTemp * xTemp;
  }
}

void func(NUMBER* x, NUMBER* y) {

  NUMBER passiveValue = codi::TypeTraits<NUMBER>::getBaseValue(x[0]);
  codi::StatementPushHelper<NUMBER> ph;

  // two valid dependencies
  ph.startPushStatement();
  ph.pushArgument(x[0], 101.0);
  ph.pushArgument(x[1], 102.0);
  ph.endPushStatement(y[0], 1.0);

  // one invalid dependencie jac == 0
  ph.startPushStatement();
  ph.pushArgument(x[0], 201.0);
  ph.pushArgument(x[1], 0.0);
  ph.endPushStatement(y[1], 2.0);

  // one invalid dependencie index == 0
  ph.startPushStatement();
  ph.pushArgument(x[0], 301.0);
  ph.pushArgument(passiveValue, 302.0);
  ph.endPushStatement(y[2], 3.0);

  // one invalid dependencie jac == inf
  ph.startPushStatement();
  ph.pushArgument(x[0], 401.0);
  ph.pushArgument(x[1], std::numeric_limits<NUMBER::Real>::infinity());
  ph.endPushStatement(y[3], 4.0);

  // one invalid dependencie jac == NaN
  ph.startPushStatement();
  ph.pushArgument(x[0], 501.0);
  ph.pushArgument(x[1], std::numeric_limits<NUMBER::Real>::quiet_NaN());
  ph.endPushStatement(y[4], 5.0);

  // two invalid dependencie jac == 0, jac == NaN
  ph.startPushStatement();
  ph.pushArgument(x[0], 0.0);
  ph.pushArgument(x[1], std::numeric_limits<NUMBER::Real>::quiet_NaN());
  ph.endPushStatement(y[5], 6.0);

  std::vector<NUMBER> inputData;
  inputData.push_back(x[0]);
  inputData.push_back(x[1]);

  std::vector<NUMBER::Real> jacData;
  jacData.push_back(701.0);
  jacData.push_back(702.0);

  // iterator push
  ph.pushStatement(y[6], 7.0, inputData.begin(), inputData.end(), jacData.begin());

  // vector push
  ph.pushStatement(y[7], 8.0, inputData, jacData, 2);
}
