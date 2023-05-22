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

#include "../tests/allTests.hpp"
#include "driverInterface.hpp"

template<typename T_Number>
struct DriverBase : public DriverInterface<T_Number> {
  public:

    using Number = CODI_DECLARE_DEFAULT(T_Number, double);

    virtual void createAllTests(TestVector<Number>& tests) = 0;

  private:

    std::string name;

  public:

    DriverBase(std::string const& name) : name(name) {}

    std::string getName() {
      return name;
    }

    TestVector<Number> getTestInfos() {
      TestVector<Number> testInfos;

      createAllTests(testInfos);

      return testInfos;
    }

  protected:
    template<typename Number>
    void createTests(TestVector<Number>& tests) {
      (void)tests;
    }

    template<typename Number, typename Test, typename... Args>
    void createTests(TestVector<Number>& tests) {
      tests.push_back(TestInfo<Number>(new Test(), Test::template func<Number>));

      createTests<Number, Args...>(tests);
    }

    void prepare(Number* x, Number* y, int curPoint, TestInterface* test, FILE* out) {
      fprintf(out, "Point %d : {", curPoint);

      for (int i = 0; i < test->getInputCount(); ++i) {
        if (i != 0) {
          fprintf(out, ", ");
        }
        double val = test->getEvalPoint(curPoint, i);
        fprintf(out, "%f", val);

        x[i] = (Number)(val);
      }
      fprintf(out, "}\n");

      for (int i = 0; i < test->getOutputCount(); ++i) {
        y[i] = 0.0;
      }
    }
};
