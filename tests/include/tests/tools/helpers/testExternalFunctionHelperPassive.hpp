/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2021 Chair for Scientific Computing (SciComp), TU Kaiserslautern
 * Homepage: http://www.scicomp.uni-kl.de
 * Contact:  Prof. Nicolas R. Gauger (codi@scicomp.uni-kl.de)
 *
 * Lead developers: Max Sagebaum, Johannes Blühdorn (SciComp, TU Kaiserslautern)
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
 *  - SciComp, TU Kaiserslautern:
 *    - Max Sagebaum
 *    - Johannes Blühdorn
 *    - Former members:
 *      - Tim Albring
 */
#pragma once

#include "../../../testInterface.hpp"
#include "multiplyExternalFunctionHelper.hpp"

struct TestExternalFunctionHelperPassive : public TestInterface {
  public:
    NAME("ExternalFunctionHelperPassive")
    IN(2)
    OUT(1)
    POINTS(1) = {{2.0, 3.0}};

    static int constexpr ITER = 5;

    template<typename Number>
    static void func(Number* x, Number* y) {
      Number w[ITER];

      w[0] = x[0];
      for (int i = 1; i < ITER; ++i) {
        w[i] = MultiplyExternalFunctionHelper<Number>::create(x[1], w[i - 1], true);
      }

      y[0] = w[ITER - 1] * w[ITER - 1];
    }
};
