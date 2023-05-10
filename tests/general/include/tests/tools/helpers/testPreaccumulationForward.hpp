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

struct TestPreaccumulationForward : public TestInterface {
  public:
    NAME("PreaccumulationForward")
    IN(2)
    OUT(4)
    POINTS(1) = {{1.0, 0.5}};

    template<typename Number>
    static void evalFunc(Number* x, Number* y) {
      y[0] = x[0];
      y[1] = x[1];
      for (int i = 0; i < 1000; ++i) {
        Number xTemp = y[0];
        Number yTemp = y[1];

        Number xSqr = xTemp * xTemp;
        Number ySqr = yTemp * yTemp;

        y[0] = xSqr - ySqr - 0.65;
        y[1] = 2.0 * yTemp * xTemp;
      }

      y[2] = x[0] * x[0];
      y[3] = x[1] * x[1];
    }

    template<typename Number>
    static void func(Number* x, Number* y) {
      codi::PreaccumulationHelper<Number> ph;

      ph.start(x[0], x[1]);

      evalFunc(x, y);

      ph.finish(false, y[0], y[1], y[2], y[3]);
    }
};
