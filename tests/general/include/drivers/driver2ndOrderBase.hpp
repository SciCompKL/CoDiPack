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

#include <codi/tools/data/hessian.hpp>

#include "../output.hpp"
#include "driverBase.hpp"

template<typename T_Number>
struct Driver2ndOrderBase : public DriverBase<T_Number> {
  public:

    using Number = CODI_DECLARE_DEFAULT(T_Number, double);
    using Base = DriverBase<Number>;

    virtual void evaluateHessian(TestInfo<Number>& info, Number* x, size_t inputs, Number* y, size_t outputs,
                                 codi::Hessian<double>& hes) = 0;

  public:

    Driver2ndOrderBase(std::string const& name) : Base(name) {}

    void runTest(TestInfo<Number>& info, FILE* out) {
      TestInterface* test = info.test;

      int evalPoints = test->getEvalPointsCount();
      int inputs = test->getInputCount();
      int outputs = test->getOutputCount();
      Number* x = new Number[inputs];
      Number* y = new Number[outputs];

      codi::Hessian<double> hes(outputs, inputs);

      for (int curPoint = 0; curPoint < evalPoints; ++curPoint) {
        Base::prepare(x, y, curPoint, test, out);

        evaluateHessian(info, x, inputs, y, outputs, hes);

        writeOutputHessian(out, hes);
      }

      delete[] x;
      delete[] y;
    }
};
