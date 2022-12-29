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

#include "../driver2ndOrderBase.hpp"

#include DRIVER_TESTS_INC

struct CoDiEvalHelper2ndOrder : public Driver2ndOrderBase<CODI_TYPE> {
  public:

    using Number = CODI_TYPE;

    using Base = Driver2ndOrderBase;

    using Gradient = Number::Gradient;

    CoDiEvalHelper2ndOrder() : Base(CODI_TO_STRING(CODI_TYPE_NAME)) {}

    void createAllTests(TestVector<Number>& tests) {
      createTests<Number, DRIVER_TESTS>(tests);
    }

    void evaluateHessian(TestInfo<Number>& info, Number* x, size_t inputs, Number* /*y*/, size_t outputs,
                         codi::Hessian<double>& hes) {
      std::vector<double> xVec(inputs);

      for (size_t i = 0; i < inputs; ++i) {
        xVec[i] = codi::RealTraits::getPassiveValue(x[i]);
      }

      auto evalFunc = [&](std::vector<Number>& x, std::vector<Number>& y) { info.func(x.data(), y.data()); };

      auto handle = codi::EvaluationHelper::template createHandle<Number>(evalFunc, outputs, inputs);

      codi::EvaluationHelper::evalHandleHessian(handle, xVec, hes);

      // evaluate a second time to force at least one tape reset.
      codi::EvaluationHelper::evalHandleHessian(handle, xVec, hes);
    }
};
