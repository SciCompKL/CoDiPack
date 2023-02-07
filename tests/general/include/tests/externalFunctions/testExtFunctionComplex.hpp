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

#include "../../testInterface.hpp"
#include "../expressions/complex/complexTestHelpers.hpp"
#include "multiplyExternalFunction.hpp"

struct TestExtFunctionComplex : public TestInterface {
  public:

    template<typename T>
    using Complex = TestComplex<T>;

    static int constexpr in_complex_count = 2;
    static int constexpr out_complex_count = 1;

    NAME("ExtFunctionComplex")
    IN(in_complex_count * 2)
    OUT(out_complex_count * 2)
    POINTS(1) =  // clang-format off
    {
      {-10.0,   5, -2.5, 1.25}
    };  // clang-format on

    template<typename Number>
    static void func(Number* x, Number* y) {
      using C = Complex<Number>;
      C xC[in_complex_count];
      C yC[out_complex_count];

      assignToComplex(xC, x, in_complex_count);

      yC[0] = MultiplyExternalFunction<C, typename Number::Tape>::create(xC[0], xC[1], Number::getTape());

      assignToReal(y, yC, out_complex_count);
    }
};
