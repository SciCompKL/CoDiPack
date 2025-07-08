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

#include "../../testInterface.hpp"
#include "../expressions/complex/complexTestHelpers.hpp"

struct TestDataExtractionTraits : public TestInterface {
  public:

    template<typename T>
    using Complex = TestComplex<T>;

    static int constexpr in_complex_count = 1;
    static int constexpr out_complex_count = 1;

    NAME("DataExtractionTraits")
    IN(in_complex_count * 2)
    OUT(out_complex_count * 2)
    POINTS(1) =  // clang-format off
    {
      {-10.0,   5}
    };  // clang-format on

    template<typename Number>
    static void func(Number* x, Number* y) {
      using C = Complex<Number>;
      C xC[in_complex_count];
      C yC[out_complex_count];

      assignToComplex(xC, x, in_complex_count);

      using DE = codi::RealTraits::DataExtraction<C>;

      // Copy over the value with the data extraction
      DE::setValue(yC[0], DE::getValue(xC[0]));
      DE::setIdentifier(yC[0], DE::getIdentifier(xC[0]));

      // Delete the old value
      DE::setValue(xC[0], typename DE::Real{});
      DE::setIdentifier(xC[0], typename DE::Identifier{});

      assignToReal(y, yC, out_complex_count);
    }
};
