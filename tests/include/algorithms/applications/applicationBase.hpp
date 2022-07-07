/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2022 Chair for Scientific Computing (SciComp), TU Kaiserslautern
 * Homepage: http://www.scicomp.uni-kl.de
 * Contact:  Prof. Nicolas R. Gauger (codi@scicomp.uni-kl.de)
 *
 * Lead developers: Max Sagebaum, Johannes Blühdorn (SciComp, TU Kaiserslautern)
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
 *  - SciComp, TU Kaiserslautern:
 *    - Max Sagebaum
 *    - Johannes Blühdorn
 *    - Former members:
 *      - Tim Albring
 */

#pragma once

#include <cstdlib>
#include <iostream>
#include <fstream>

#include <codi.hpp>

#include "applicationBaseSettings.hpp"
#include "checkpointManager.hpp"
#include "io.hpp"

template<typename T_Type, typename T_Impl>
struct TestApplicationBase : public codi::algorithms::ApplicationInterface<T_Type> {
  public:
    using Type = CODI_DD(T_Type, CODI_T(codi::LhsExpressionInterface<double, double, CODI_ANY, CODI_ANY>));
    using Impl = CODI_DD(T_Impl, CODI_T(codi::algorithms::ApplicationInterface<CODI_ANY>));

    using Real = typename Type::Real;

  private:
    int curIteration;

    std::ofstream out;

    TestCheckpointManager<Type, Impl> testCPM;
    TestIO<Type, Impl> testIO;

  public:

    TestApplicationBaseSettings generalSettings;

    TestApplicationBase(std::string file, Impl* impl) :
      curIteration(0),
      out(file),
      testCPM(impl),
      testIO(impl),
      generalSettings()
    {}

    codi::algorithms::Residuum<Real> residuumY(std::vector<Real> const& v1, std::vector<Real> const& v2) {
      codi::algorithms::Residuum<Real> res{};
      res.lMax = -1e300;

      for(size_t i = 0; i < v1.size(); i += 1) {
        Real diff = abs(v1[i] - v2[i]);
        res.l1 += diff;
        res.l2 += diff * diff;
        if(res.lMax < diff) {
          res.lMax = diff;
          res.lMaxPos = i;
        }
      }

      res.l2 = sqrt(res.l2);

      return res;
    }

    codi::algorithms::Residuum<Real> residuumX(std::vector<Real> const& v1, std::vector<Real> const& v2) {
      return residuumY(v1, v2);
    }

    codi::algorithms::Residuum<Real> residuumP(std::vector<Real> const& v1, std::vector<Real> const& v2) {
      return residuumY(v1, v2);
    }

    TestCheckpointManager<Type, Impl>* getCheckpointInterface() {
      return &testCPM;
    }

    TestIO<Type, Impl>* getIOInterface() {
      return &testIO;
    }

    codi::algorithms::ApplicationHints getHints() {
      return codi::algorithms::ApplicationHints::NONE();
    }

    int getIteration() {
      return curIteration;
    }

    void print(std::string const& line) {
      if(generalSettings.writeToStdout) {
        std::cout << line;
      }
      out << line;
    }

    bool isStop() {
      return false;
    }

    void setIteration(int iteration) {
      curIteration = iteration;
    }

    void setOutputFile(std::string const& file) {
      out.close();
      out.open(file);
    }
};
