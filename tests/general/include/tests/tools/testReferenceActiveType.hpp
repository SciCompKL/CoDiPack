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
#include <codi.hpp>
#include <type_traits>
#include <vector>

#include "../../testInterface.hpp"

struct TestReferenceActiveType : public TestInterface {
  public:
    NAME("ReferenceActiveType")
    IN(1)
    OUT(1)
    POINTS(1) = {// clang-format off
      {0.5}
    };  // clang-format on

#if !defined(DOUBLE)
    template<typename Number>
    using RefReal = codi::ReferenceActiveType<Number>;
#else
    template<typename Number>
    using RefReal = Number;
#endif

    template<typename Number>
    static void func(Number* x, Number* y) {
      RefReal<Number> xRef = x[0];
      y[0] = 3.0 * xRef * xRef * xRef * xRef + 5.0 * xRef * xRef * xRef - 3.0 * xRef * xRef + 2.0 * xRef - 4.0;
    }
};
