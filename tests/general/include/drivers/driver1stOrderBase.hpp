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

#include <codi/tools/data/jacobian.hpp>

#include "../output.hpp"
#include "driverBase.hpp"

template<typename T_Number>
struct Driver1stOrderBase : public DriverBase<T_Number> {
  public:

    using Number = CODI_DECLARE_DEFAULT(T_Number, double);
    using Base = DriverBase<Number>;

    virtual void evaluateJacobian(TestInfo<Number>& info, Number* x, size_t inputs, Number* y, size_t outputs,
                                  codi::Jacobian<double>& jac) = 0;

  public:

    Driver1stOrderBase(std::string const& name) : Base(name) {}

    void runTest(TestInfo<Number>& info, FILE* out) {
      TestInterface* test = info.test;

      int evalPoints = test->getEvalPointsCount();
      int inputs = test->getInputCount();
      int outputs = test->getOutputCount();
      Number* x = new Number[inputs];
      Number* y = new Number[outputs];

      codi::Jacobian<double> jac(outputs, inputs);

      for (int curPoint = 0; curPoint < evalPoints; ++curPoint) {
        Base::prepare(x, y, curPoint, test, out);

        evaluateJacobian(info, x, inputs, y, outputs, jac);

        writeOutputJacobian(out, jac);
      }

      delete[] x;
      delete[] y;
    }
};
