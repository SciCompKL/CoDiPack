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

#include <codi.hpp>

#include "../../../testInterface.hpp"
#include "complexTestHelpers.hpp"

struct TestComplexTwoArgumentExpr1 : public TestInterface {
  public:

    template<typename T>
    using Complex = TestComplex<T>;

    static int constexpr in_complex_count = 2;
    static int constexpr out_complex_count = 48;
    static int constexpr out_real_count = 0;

    static int constexpr out_real_offset = out_complex_count * 2;

    NAME("ComplexTwoArgumentExpr1")
    IN(in_complex_count * 2)
    OUT(out_complex_count * 2 + out_real_count)
    POINTS(81) =  // clang-format off
    {
      {-10.0,  -5.0,  -2.0,  -1.0},
      {-10.0,  -5.0,  -2.0,   0.0},
      {-10.0,  -5.0,  -2.0,   1.0},
      {-10.0,  -5.0,   0.0,  -1.0},
      {-10.0,  -5.0,   0.0,   0.0},
      {-10.0,  -5.0,   0.0,   1.0},
      {-10.0,  -5.0,   2.0,  -1.0},
      {-10.0,  -5.0,   2.0,   0.0},
      {-10.0,  -5.0,   2.0,   1.0},
      {-10.0,   0.0,  -2.0,  -1.0},
      {-10.0,   0.0,  -2.0,   0.0},
      {-10.0,   0.0,  -2.0,   1.0},
      {-10.0,   0.0,   0.0,  -1.0},
      {-10.0,   0.0,   0.0,   0.0},
      {-10.0,   0.0,   0.0,   1.0},
      {-10.0,   0.0,   2.0,  -1.0},
      {-10.0,   0.0,   2.0,   0.0},
      {-10.0,   0.0,   2.0,   1.0},
      {-10.0,   5.0,  -2.0,  -1.0},
      {-10.0,   5.0,  -2.0,   0.0},
      {-10.0,   5.0,  -2.0,   1.0},
      {-10.0,   5.0,   0.0,  -1.0},
      {-10.0,   5.0,   0.0,   0.0},
      {-10.0,   5.0,   0.0,   1.0},
      {-10.0,   5.0,   2.0,  -1.0},
      {-10.0,   5.0,   2.0,   0.0},
      {-10.0,   5.0,   2.0,   1.0},
      {  0.0,  -5.0,  -2.0,  -1.0},
      {  0.0,  -5.0,  -2.0,   0.0},
      {  0.0,  -5.0,  -2.0,   1.0},
      {  0.0,  -5.0,   0.0,  -1.0},
      {  0.0,  -5.0,   0.0,   0.0},
      {  0.0,  -5.0,   0.0,   1.0},
      {  0.0,  -5.0,   2.0,  -1.0},
      {  0.0,  -5.0,   2.0,   0.0},
      {  0.0,  -5.0,   2.0,   1.0},
      {  0.0,   0.0,  -2.0,  -1.0},
      {  0.0,   0.0,  -2.0,   0.0},
      {  0.0,   0.0,  -2.0,   1.0},
      {  0.0,   0.0,   0.0,  -1.0},
      {  0.0,   0.0,   0.0,   0.0},
      {  0.0,   0.0,   0.0,   1.0},
      {  0.0,   0.0,   2.0,  -1.0},
      {  0.0,   0.0,   2.0,   0.0},
      {  0.0,   0.0,   2.0,   1.0},
      {  0.0,   5.0,  -2.0,  -1.0},
      {  0.0,   5.0,  -2.0,   0.0},
      {  0.0,   5.0,  -2.0,   1.0},
      {  0.0,   5.0,   0.0,  -1.0},
      {  0.0,   5.0,   0.0,   0.0},
      {  0.0,   5.0,   0.0,   1.0},
      {  0.0,   5.0,   2.0,  -1.0},
      {  0.0,   5.0,   2.0,   0.0},
      {  0.0,   5.0,   2.0,   1.0},
      { 10.0,  -5.0,  -2.0,  -1.0},
      { 10.0,  -5.0,  -2.0,   0.0},
      { 10.0,  -5.0,  -2.0,   1.0},
      { 10.0,  -5.0,   0.0,  -1.0},
      { 10.0,  -5.0,   0.0,   0.0},
      { 10.0,  -5.0,   0.0,   1.0},
      { 10.0,  -5.0,   2.0,  -1.0},
      { 10.0,  -5.0,   2.0,   0.0},
      { 10.0,  -5.0,   2.0,   1.0},
      { 10.0,   0.0,  -2.0,  -1.0},
      { 10.0,   0.0,  -2.0,   0.0},
      { 10.0,   0.0,  -2.0,   1.0},
      { 10.0,   0.0,   0.0,  -1.0},
      { 10.0,   0.0,   0.0,   0.0},
      { 10.0,   0.0,   0.0,   1.0},
      { 10.0,   0.0,   2.0,  -1.0},
      { 10.0,   0.0,   2.0,   0.0},
      { 10.0,   0.0,   2.0,   1.0},
      { 10.0,   5.0,  -2.0,  -1.0},
      { 10.0,   5.0,  -2.0,   0.0},
      { 10.0,   5.0,  -2.0,   1.0},
      { 10.0,   5.0,   0.0,  -1.0},
      { 10.0,   5.0,   0.0,   0.0},
      { 10.0,   5.0,   0.0,   1.0},
      { 10.0,   5.0,   2.0,  -1.0},
      { 10.0,   5.0,   2.0,   0.0},
      { 10.0,   5.0,   2.0,   1.0}
    };  // clang-format on

    template<typename Number>
    static void func(Number* x, Number* y) {
      using C = Complex<Number>;
      C xC[in_complex_count];
      C yC[out_complex_count];

      assignToComplex(xC, x, in_complex_count);

      // clang-format off
       yC[0] =         xC[0]  +         xC[1];
       yC[1] =         xC[0]  + passive(xC[1]);
       yC[2] = passive(xC[0]) +         xC[1];
       yC[3] =         xC[0]  +          x[2];
       yC[4] =         xC[0]  +  passive(x[2]);
       yC[5] = passive(xC[0]) +          x[2];
       yC[6] =          x[0]  +         xC[1];
       yC[7] =          x[0]  + passive(xC[1]);
       yC[8] =  passive(x[0]) +         xC[1];
       yC[9] =         xC[0]  -         xC[1];
      yC[10] =         xC[0]  - passive(xC[1]);
      yC[11] = passive(xC[0]) -         xC[1];
      yC[12] =         xC[0]  -          x[2];
      yC[13] =         xC[0]  -  passive(x[2]);
      yC[14] = passive(xC[0]) -          x[2];
      yC[15] =          x[0]  -         xC[1];
      yC[16] =          x[0]  - passive(xC[1]);
      yC[17] =  passive(x[0]) -         xC[1];
      yC[18] =         xC[0]  *         xC[1];
      yC[19] =         xC[0]  * passive(xC[1]);
      yC[20] = passive(xC[0]) *         xC[1];
      yC[21] =         xC[0]  *          x[2];
      yC[22] =         xC[0]  *  passive(x[2]);
      yC[23] = passive(xC[0]) *          x[2];
      yC[24] =          x[0]  *         xC[1];
      yC[25] =          x[0]  * passive(xC[1]);
      yC[26] =  passive(x[0]) *         xC[1];
      yC[27] =         xC[0]  /         xC[1];
      yC[28] =         xC[0]  / passive(xC[1]);
      yC[29] = passive(xC[0]) /         xC[1];
      yC[30] =         xC[0]  /          x[2];
      yC[31] =         xC[0]  /  passive(x[2]);
      yC[32] = passive(xC[0]) /          x[2];
      yC[33] =          x[0]  /         xC[1];
      yC[34] =          x[0]  / passive(xC[1]);
      yC[35] =  passive(x[0]) /         xC[1];
      if(0.0 != xC[0]) {
        yC[36] = pow(        xC[0] ,         xC[1]);
        yC[37] = pow(        xC[0] , passive(xC[1]));
        yC[38] = pow(passive(xC[0]),         xC[1]);
        yC[39] = pow(        xC[0] ,          x[2]);
        yC[40] = pow(        xC[0] ,  passive(x[2]));
        yC[41] = pow(passive(xC[0]),          x[2]);
        yC[42] = pow(         x[0] ,         xC[1]);
        yC[43] = pow(         x[0] , passive(xC[1]));
        yC[44] = pow( passive(x[0]),         xC[1]);
      } else {
        yC[36] = 0.0;
        yC[37] = 0.0;
        yC[38] = 0.0;
        yC[39] = 0.0;
        yC[40] = 0.0;
        yC[41] = 0.0;
        yC[42] = 0.0;
        yC[43] = 0.0;
        yC[44] = 0.0;
      }
      if(x[0] >= 0.0) {
        yC[45] = polar(        x[0],         x[1]);
        yC[46] = polar(passive(x[0]),        x[1]);
        yC[47] = polar(        x[0], passive(x[1]));
      } else {
        yC[45] = 0.0;
        yC[46] = 0.0;
        yC[47] = 0.0;
      }
      // clang-format on

      assignToReal(y, yC, out_complex_count);
    }
};
