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
#include "../../../testInterface.hpp"

struct TestPreaccumulationForwardInvalidAdjoint : public TestInterface {
  public:
    NAME("PreaccumulationForwardInvalidAdjoint")
    IN(2)
    OUT(4)
    POINTS(1) = {{1.0, 0.5}};

    template<typename Number>
    static void evalFunc(Number* x, Number* y) {
      Number temp1 = x[0] * x[1];
      Number temp2 = x[0] / x[1];
      Number temp3 = x[0] + x[1];
      Number temp4 = x[0] - x[1];
      Number temp5 = temp1 * temp3;
      Number temp6 = temp2 * temp4;

      y[0] = temp5 * temp5;
      y[1] = temp6 * temp6;
      y[2] = temp5 * temp5;
      y[3] = temp6 * temp6;
    }

    template<typename Number>
    static void func(Number* x, Number* y) {
      codi::PreaccumulationHelper<Number> ph;

      ph.start(x[0], x[1]);

      evalFunc(x, y);

      ph.finish(false, y[0], y[1], y[2], y[3]);

      Number temp1 = y[0] + y[1];
      Number temp2 = y[2] + y[3];

      y[0] = temp1 + temp2;
      y[1] = temp1 - temp2;
      y[2] = temp1 * temp2;
      y[3] = temp1 / temp2;
    }
};
