/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2026 Chair for Scientific Computing (SciComp), RPTU University Kaiserslautern-Landau
 * Homepage: http://scicomp.rptu.de
 * Contact:  Prof. Nicolas R. Gauger (codi@scicomp.uni-kl.de)
 *
 * Lead developers: Max Sagebaum, Johannes Blühdorn (SciComp, RPTU University Kaiserslautern-Landau)
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
 *  - SciComp, RPTU University Kaiserslautern-Landau:
 *    - Max Sagebaum
 *    - Johannes Blühdorn
 *    - Former members:
 *      - Tim Albring
 */

#pragma once

#include "../../../testInterface.hpp"
#include "complexTestHelpers.hpp"

struct TestComplexOneArgumentExpr1 : public TestInterface {
  public:

    template<typename T>
    using Complex = TestComplex<T>;

    static int constexpr in_complex_count = 1;
    static int constexpr out_complex_count = 16;
    static int constexpr out_real_count = 7;

    static int constexpr out_real_offset = out_complex_count * 2;

    NAME("ComplexOneArgumentExpr1")
    IN(in_complex_count * 2)
    OUT(out_complex_count * 2 + out_real_count)
    POINTS(25) =  // clang-format off
    {
      {-10.0,   -10},
      {-10.0,    -5},
      {-10.0,     0},
      {-10.0,     5},
      {-10.0,    10},
      { -5.0,   -10},
      { -5.0,    -5},
      { -5.0,     0},
      { -5.0,     5},
      { -5.0,    10},
      {  0.0,   -10},
      {  0.0,    -5},
      {  0.0,     0},
      {  0.0,     5},
      {  0.0,    10},
      {  5.0,   -10},
      {  5.0,    -5},
      {  5.0,     0},
      {  5.0,     5},
      {  5.0,    10},
      { 10.0,   -10},
      { 10.0,    -5},
      { 10.0,     0},
      { 10.0,     5},
      { 10.0,    10}
    };  // clang-format on

    template<typename Number>
    static void func(Number* x, Number* y) {
      int constexpr realOffset = out_real_offset;
      using C = Complex<Number>;
      C xC[in_complex_count];
      C yC[out_complex_count];

      assignToComplex(xC, x, in_complex_count);

      yC[0] = conj(xC[0]);               // R x R
      yC[1] = proj(xC[0]);               // R x R
      yC[2] = exp(xC[0]);                // R x R
      yC[3] = log(xC[0]);                // R x R \ {0, 0}
      yC[4] = log10(xC[0]);              // R x R \ {0, 0}
      yC[5] = sin(xC[0]);                // R x R
      yC[6] = cos(xC[0]);                // R x R
      yC[7] = tan(xC[0]);                // R x R \ {{(1/2 + i) * PI, 0} | for all i in Z}
      yC[8] = atan(xC[0]);               // R x R \ {{0, 1}, {0, -1}
      yC[9] = sinh(xC[0]);               // R x R
      yC[10] = cosh(xC[0]);              // R x R
      yC[11] = tanh(xC[0]);              // R x R
      yC[12] = acosh(xC[0]);             // R x R \ {{1, 0}, {-1, 0}}
      yC[13] = atanh(xC[0]);             // R x R \ {{1, 0}, {-1, 0}}
      yC[14] = +xC[0];                   // R x R
      yC[15] = -xC[0];                   // R x R
      y[realOffset + 0] = xC[0].real();  // R x R
      y[realOffset + 1] = xC[0].imag();  // R x R
      y[realOffset + 2] = real(xC[0]);   // R x R
      y[realOffset + 3] = imag(xC[0]);   // R x R
      y[realOffset + 4] = abs(xC[0]);    // R x R \ {0,0}
      y[realOffset + 5] = arg(xC[0]);    // R x R \ {0,0}
      y[realOffset + 6] = norm(xC[0]);   // R x R

      assignToReal(y, yC, out_complex_count);
    }
};
