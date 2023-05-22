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
#include <codi/tools/data/jacobian.hpp>

#include "../driver1stOrderBase.hpp"

#include DRIVER_TESTS_INC

struct CoDiForward1stOrder : public Driver1stOrderBase<CODI_TYPE> {
  public:

    using Number = CODI_TYPE;

    using Base = Driver1stOrderBase<Number>;

    CoDiForward1stOrder() : Base(CODI_TO_STRING(CODI_TYPE_NAME)) {}

    void createAllTests(TestVector<Number>& tests) {
      createTests<Number, DRIVER_TESTS>(tests);
    }

    void evaluateJacobian(TestInfo<Number>& info, Number* x, size_t inputs, Number* y, size_t outputs,
                          codi::Jacobian<double>& jac) {
      using Gradient = typename Number::Gradient;
      size_t constexpr gradDim = codi::GradientTraits::dim<Gradient>();

      size_t runs = inputs / gradDim;
      if (inputs % gradDim != 0) {
        runs += 1;
      }
      for (size_t curIn = 0; curIn < runs; ++curIn) {
        size_t curSize = gradDim;
        if ((curIn + 1) * gradDim > (size_t)inputs) {
          curSize = inputs % gradDim;
        }
        for (size_t curDim = 0; curDim < curSize; ++curDim) {
          codi::GradientTraits::at(x[curIn * gradDim + curDim].gradient(), curDim) = 1.0;
        }

        for (size_t i = 0; i < outputs; ++i) {
          y[i].setGradient(Gradient());
        }

        info.func(x, y);

        for (size_t curDim = 0; curDim < curSize; ++curDim) {
          for (size_t curOut = 0; curOut < outputs; ++curOut) {
#ifdef SECOND_ORDER
            jac(curOut, curIn * gradDim + curDim) = codi::GradientTraits::at(y[curOut].getGradient(), curDim).value();
#else
            jac(curOut, curIn * gradDim + curDim) = codi::GradientTraits::at(y[curOut].getGradient(), curDim);
#endif
          }
        }

        for (size_t curDim = 0; curDim < curSize; ++curDim) {
          x[curIn * gradDim + curDim].setGradient(Gradient());
        }
      }
    }
};
