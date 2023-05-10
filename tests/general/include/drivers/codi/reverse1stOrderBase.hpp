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

struct CoDiReverse1stOrderBase : public Driver1stOrderBase<CODI_TYPE> {
  public:

    using Number = CODI_DECLARE_DEFAULT(
        CODI_TYPE, CODI_TEMPLATE(codi::LhsExpressionInterface<double, double, CODI_ANY, CODI_ANY>));

    using Tape = CODI_DD(typename Number::Tape, CODI_T(codi::FullTapeInterface<double, double, int, CODI_ANY>));
    using Base = Driver1stOrderBase<Number>;

    using Gradient = Number::Gradient;

    CoDiReverse1stOrderBase() : Base(CODI_TO_STRING(CODI_TYPE_NAME)) {}

    virtual Gradient& accessGradient(Number& value) = 0;
    virtual void cleanup() = 0;
    virtual void evaluate() = 0;
    virtual void prepare() = 0;

    void createAllTests(TestVector<Number>& tests) {
      createTests<Number, DRIVER_TESTS>(tests);
    }

    void evaluateJacobian(TestInfo<Number>& info, Number* x, size_t inputs, Number* y, size_t outputs,
                          codi::Jacobian<double>& jac) {
      using Gradient = typename Number::Gradient;
      size_t constexpr gradDim = codi::GradientTraits::dim<Gradient>();

      Tape& tape = Number::getTape();

      setTapeSizes(tape);

      size_t runs = outputs / gradDim;
      if (outputs % gradDim != 0) {
        runs += 1;
      }

      for (size_t curOut = 0; curOut < runs; ++curOut) {
        size_t curSize = gradDim;
        if ((curOut + 1) * gradDim > (size_t)outputs) {
          curSize = outputs % gradDim;
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

        for (size_t curDim = 0; curDim < curSize; ++curDim) {
          if (tape.isIdentifierActive(y[curOut * gradDim + curDim].getIdentifier())) {
            codi::GradientTraits::at(accessGradient(y[curOut * gradDim + curDim]), curDim) = 1.0;
          }
        }

        evaluate();

        for (size_t curDim = 0; curDim < curSize; ++curDim) {
          for (size_t curIn = 0; curIn < inputs; ++curIn) {
#ifdef SECOND_ORDER
            jac(curOut * gradDim + curDim, curIn) = codi::GradientTraits::at(accessGradient(x[curIn]), curDim).value();
#else
            jac(curOut * gradDim + curDim, curIn) = codi::GradientTraits::at(accessGradient(x[curIn]), curDim);
#endif
          }
        }

        tape.reset();

        cleanup();
      }
    }

  protected:

    void setTapeSizes(Tape& tape) {
      // Set sizes for Jacobian tapes
      if (tape.hasParameter(codi::TapeParameters::JacobianSize)) {
        tape.setParameter(codi::TapeParameters::JacobianSize, 10000);
      }
      if (tape.hasParameter(codi::TapeParameters::StatementSize)) {
        tape.setParameter(codi::TapeParameters::StatementSize, 10000);
      }
      if (tape.hasParameter(codi::TapeParameters::ExternalFunctionsSize)) {
        tape.setParameter(codi::TapeParameters::ExternalFunctionsSize, 10000);
      }
    }
};
