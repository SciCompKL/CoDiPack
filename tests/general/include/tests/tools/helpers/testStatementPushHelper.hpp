/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2023 Chair for Scientific Computing (SciComp), University of Kaiserslautern-Landau
 * Homepage: http://www.scicomp.uni-kl.de
 * Contact:  Prof. Nicolas R. Gauger (codi@scicomp.uni-kl.de)
 *
 * Lead developers: Max Sagebaum, Johannes Blühdorn (SciComp, University of Kaiserslautern-Landau)
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
 * For other licensing options please contact us.
 *
 * Authors:
 *  - SciComp, University of Kaiserslautern-Landau:
 *    - Max Sagebaum
 *    - Johannes Blühdorn
 *    - Former members:
 *      - Tim Albring
 */
#include <codi.hpp>
#include <type_traits>
#include <vector>

#include "../../../testInterface.hpp"

struct TestStatementPushHelper : public TestInterface {
  public:
    NAME("StatementPushHelper")
    IN(2)
    OUT(8)
    POINTS(1) = {// clang-format off
      {  1.0,     0.5}
    };  // clang-format on

    template<typename Number>
    static void func(Number* x, Number* y) {
      using Real = codi::RealTraits::PassiveReal<Number>;
      Number passiveValue = codi::RealTraits::getPassiveValue(x[0]);

      codi::StatementPushHelper<Number> ph;

      // two valid dependencies
      ph.startPushStatement();
      ph.pushArgument(x[0], 101.0);
      ph.pushArgument(x[1], 102.0);
      ph.endPushStatement(y[0], 1.0);

      // one invalid dependency jac == 0
      ph.startPushStatement();
      ph.pushArgument(x[0], 201.0);
      ph.pushArgument(x[1], 0.0);
      ph.endPushStatement(y[1], 2.0);

      // one invalid dependency index == 0
      ph.startPushStatement();
      ph.pushArgument(x[0], 301.0);
      ph.pushArgument(passiveValue, 302.0);
      ph.endPushStatement(y[2], 3.0);

      // one invalid dependency jac == inf
      ph.startPushStatement();
      ph.pushArgument(x[0], 401.0);
      ph.pushArgument(x[1], std::numeric_limits<Real>::infinity());
      ph.endPushStatement(y[3], 4.0);

      // one invalid dependency jac == NaN
      ph.startPushStatement();
      ph.pushArgument(x[0], 501.0);
      ph.pushArgument(x[1], std::numeric_limits<Real>::quiet_NaN());
      ph.endPushStatement(y[4], 5.0);

      // two invalid dependencies jac == 0, jac == NaN
      ph.startPushStatement();
      ph.pushArgument(x[0], 0.0);
      ph.pushArgument(x[1], std::numeric_limits<Real>::quiet_NaN());
      ph.endPushStatement(y[5], 6.0);

      std::vector<Number> inputData;
      inputData.push_back(x[0]);
      inputData.push_back(x[1]);

      std::vector<Real> jacData;
      jacData.push_back(701.0);
      jacData.push_back(702.0);

      // iterator push
      ph.pushStatement(y[6], 7.0, inputData.begin(), inputData.end(), jacData.begin());

      // vector push
      ph.pushStatement(y[7], 8.0, inputData, jacData, 2);
    }
};
