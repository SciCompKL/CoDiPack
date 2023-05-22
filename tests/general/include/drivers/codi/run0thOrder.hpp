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
#pragma once

#include <codi.hpp>

#include "../driver0thOrderBase.hpp"

#include DRIVER_TESTS_INC

struct CoDi0thOrder : public Driver0thOrderBase<CODI_TYPE> {
  public:

    using Number = CODI_TYPE;

    using Base = Driver0thOrderBase<Number>;

    CoDi0thOrder() : Base(CODI_TO_STRING(CODI_TYPE_NAME)) {}

    void createAllTests(TestVector<Number>& tests) {
      createTests<Number, DRIVER_TESTS>(tests);
    }

    void evaluatePrimal(TestInfo<Number>& info, Number* x, size_t inputs, Number* y, size_t outputs,
                        std::vector<double>& primals) {
      codi::CODI_UNUSED(inputs);

      info.func(x, y);

      for (size_t curOut = 0; curOut < outputs; curOut += 1) {
        primals[curOut] = codi::RealTraits::getPassiveValue(y[curOut]);
      }
    }
};
