/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2025 Chair for Scientific Computing (SciComp), University of Kaiserslautern-Landau
 * Homepage: http://scicomp.rptu.de
 * Contact:  Prof. Nicolas R. Gauger (codi@scicomp.uni-kl.de)
 *
 * Lead developers: Max Sagebaum, Johannes Blühdorn (SciComp, University of Kaiserslautern-Landau)
 *
 * This file is part of CoDiPack (http://scicomp.rptu.de/software/codi).
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

#pragma once

#include "../../../testInterface.hpp"

template<typename Impl>
struct BasePreaccumulationPassiveValue : public TestInterface {
  public:
    IN(2)
    OUT(2)
    POINTS(1) = { {1.0, 0.5} };

    template<typename Number>
    static void evalFunc(Number* x, Number* y) {
      y[0] = codi::RealTraits::getPassiveValue(x[0]);  // kill x dependency
      y[1] = x[1];

      Number two = 2.0;
      Number zeroSix = 0.65;
      for (int i = 0; i < 5; ++i) {
        Number xTemp = y[0];
        Number yTemp = y[1];

        y[0] = xTemp * xTemp - yTemp * yTemp - zeroSix;
        y[1] = two * yTemp * xTemp;
      }
    }

    template<typename Number>
    static void func(Number* x, Number* y) {
      codi::PreaccumulationHelper<Number> ph;

      ph.start(x[0], x[1]);

      evalFunc(x, y);

      Impl::finish(ph, y);
    }
};
