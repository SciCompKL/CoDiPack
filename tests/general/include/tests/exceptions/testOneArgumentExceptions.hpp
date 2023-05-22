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
#define _USE_MATH_DEFINES
#include <cmath>

#include "../../testInterface.hpp"

struct TestOneArgumentExceptions : public TestInterface {
  public:
    NAME("OneArgumentExceptions")
    IN(1)
    OUT(8)
    POINTS(6) =  // clang-format off
    {
      {0.5 * M_PI},
      { 0.000},
      { 1.000},
      {-1.000},
      { 2.000},
      {-2.000},
    };  // clang-format on

    template<typename Number>
    static void func(Number* x, Number* y) {
      y[0] = tan(x[0]);    // R \ {(1/2 + i) * PI | for all i in Z}
      y[1] = log(x[0]);    // (0, inf)
      y[2] = log10(x[0]);  // (0, inf)
      y[3] = sqrt(x[0]);   // [0 , inf)
      y[4] = atanh(x[0]);  // (-1, 1)
      y[5] = asin(x[0]);   // [-1, 1]
      y[6] = acos(x[0]);   // [-1, 1]
      y[7] = cbrt(x[0]);   // R \ {0}
    }
};
