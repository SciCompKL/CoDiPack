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

#include "reverse1stOrderBase.hpp"

#include DRIVER_TESTS_INC

struct CoDiForwardTape1stOrder : public CoDiReverse1stOrderBase {
  public:

    using Number = CODI_DECLARE_DEFAULT(
        CODI_TYPE, CODI_TEMPLATE(codi::LhsExpressionInterface<double, double, CODI_ANY, CODI_ANY>));

    using Tape = CODI_DD(typename Number::Tape, CODI_T(codi::FullTapeInterface<double, double, int, CODI_ANY>));
    using Base = CoDiReverse1stOrderBase;

    using Gradient = Number::Gradient;

    using Base::Base;

    Gradient& accessGradient(Number& value) {
      return value.gradient();
    }

    void cleanup() {}

    void evaluate() {
      Number::getTape().evaluateForward();
    }

    void prepare() {}

    void evaluateJacobian(TestInfo<Number>& info, Number* x, size_t inputs, Number* y, size_t outputs,
                          codi::Jacobian<double>& jac) {
      size_t constexpr gradDim = codi::GradientTraits::dim<Gradient>();

      Tape& tape = Number::getTape();

      Base::setTapeSizes(tape);

      size_t runs = inputs / gradDim;
      if (inputs % gradDim != 0) {
        runs += 1;
      }

      for (size_t curRun = 0; curRun < runs; ++curRun) {
        size_t curSize = gradDim;
        if ((curRun + 1) * gradDim > (size_t)inputs) {
          curSize = inputs % gradDim;
        }

        prepare();

        tape.setActive();

        for (size_t i = 0; i < inputs; ++i) {
          tape.registerInput(x[i]);
        }

        info.func(x, y);

        for (size_t i = 0; i < outputs; ++i) {
          tape.registerOutput(y[i]);
        }

        for (size_t i = 0; i < inputs; ++i) {
          tape.setPrimal(x[i].getIdentifier(), x[i].getValue());
        }

        for (size_t curDim = 0; curDim < curSize; ++curDim) {
          codi::GradientTraits::at(accessGradient(x[curRun * gradDim + curDim]), curDim) = 1.0;
        }

        evaluate();  // Forward evaluate in this case

        for (size_t curOut = 0; curOut < outputs; ++curOut) {
          for (size_t curDim = 0; curDim < curSize; ++curDim) {
#ifdef SECOND_ORDER
            jac(curOut, curRun * gradDim + curDim) =
                codi::GradientTraits::at(accessGradient(y[curOut]), curDim).value();
#else
            jac(curOut, curRun * gradDim + curDim) = codi::GradientTraits::at(accessGradient(y[curOut]), curDim);
#endif
          }
        }

        tape.reset();

        cleanup();
      }
    }
};
