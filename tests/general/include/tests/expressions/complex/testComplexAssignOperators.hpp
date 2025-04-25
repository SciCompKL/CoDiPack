/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2022 Chair for Scientific Computing (SciComp), TU Kaiserslautern
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

#include <codi.hpp>

#include "../../../testInterface.hpp"
#include "complexTestHelpers.hpp"

struct TestComplexAssignOperators : public TestInterface {
  public:

    template<typename T>
    using Complex = TestComplex<T>;

    static int constexpr in_complex_count = 2;
    static int constexpr out_complex_count = 16;
    static int constexpr out_real_count = 0;

    static int constexpr out_real_offset = out_complex_count * 2;

    NAME("ComplexAssignOperators")
    IN(in_complex_count * 2)
    OUT(out_complex_count * 2 + out_real_count)
    POINTS(5) =  // clang-format off
    {
      {10.0,  5.0,  2.0,  1.0},
      {10.0,  5.0,  2.0,  0.0},
      {10.0,  5.0,  0.0,  1.0},
      {10.0,  0.0,  2.0,  1.0},
      { 0.0,  5.0,  2.0,  1.0},
    };  // clang-format on

    template<typename Number>
    static void func(Number* x, Number* y) {
      using C = Complex<Number>;
      C xC[in_complex_count];
      C yC[out_complex_count];

      assignToComplex(xC, x, in_complex_count);

      // clang-format off
       yC[0] = xC[0];  yC[0] +=         xC[1];
       yC[1] = xC[0];  yC[1] += passive(xC[1]);
       yC[2] = xC[0];  yC[2] +=          x[2];
       yC[3] = xC[0];  yC[3] +=  passive(x[2]);
       yC[4] = xC[0];  yC[4] -=         xC[1];
       yC[5] = xC[0];  yC[5] -= passive(xC[1]);
       yC[6] = xC[0];  yC[6] -=          x[2];
       yC[7] = xC[0];  yC[7] -=  passive(x[2]);
       yC[8] = xC[0];  yC[8] *=         xC[1];
       yC[9] = xC[0];  yC[9] *= passive(xC[1]);
      yC[10] = xC[0]; yC[10] *=          x[2];
      yC[11] = xC[0]; yC[11] *=  passive(x[2]);
      yC[12] = xC[0]; yC[12] /=         xC[1];
      yC[13] = xC[0]; yC[13] /= passive(xC[1]);
      yC[14] = xC[0]; yC[14] /=          x[2];
      yC[15] = xC[0]; yC[15] /=  passive(x[2]);
      // clang-format on

      assignToReal(y, yC, out_complex_count);
    }
};
