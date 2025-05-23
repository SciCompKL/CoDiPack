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
#include "../../testInterface.hpp"

struct TestTwoArgumentExceptions : public TestInterface {
  public:
    NAME("TwoArgumentExceptions")
    IN(2)
    OUT(12)
    POINTS(4) =  // clang-format off
    {
      {0.0,   0.0},
      {1.0,   0.0},
      {0.0,  -5.0},
      {-1.0,  0.5}
    };  // clang-format on

    template<typename Number>
    static void func(Number* x, Number* y) {
      y[0] = x[0] / x[1];            // R x (R \ {0})
      y[1] = 5.00 / x[1];            // R x (R \ {0})
      y[2] = x[0] / 5.00;            // R x (R \ {0})
      y[3] = atan2(x[0], x[1]);      // R x R \ {0, 0}
      y[4] = atan2(5.00, x[1]);      // R x R \ {0, 0}
      y[5] = atan2(x[0], 5.00);      // R x R \ {0, 0}
      y[6] = remainder(x[0], x[1]);  // R x (R \ {0})
      y[7] = remainder(5.00, x[1]);  // R x (R \ {0})
      y[8] = remainder(x[0], 5.00);  // R x (R \ {0})
      y[9] = hypot(x[0], x[1]);      // R x R \ {0, 0})
      y[10] = hypot(5.0, x[1]);      // R x R \ {0, 0})
      y[11] = hypot(x[0], 5.0);      // R x R \ {0, 0})
    }
};
