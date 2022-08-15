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

#include <codi.hpp>
#include <cstdlib>
#include <fstream>
#include <iostream>

template<typename T_Type, typename T_Impl>
struct TestApplicationBase : public codi::algorithms::DefaultApplication<T_Type, T_Impl> {
  public:
    using Type = CODI_DD(T_Type, CODI_T(codi::LhsExpressionInterface<double, double, CODI_ANY, CODI_ANY>));
    using Impl = CODI_DD(T_Impl, CODI_T(codi::algorithms::ApplicationInterface<CODI_ANY>));

    using Base = codi::algorithms::DefaultApplication<Type, Impl>;

    using Real = typename Type::Real;
    using Res = typename Base::Res;

  private:
    std::ofstream out;

  public:

    bool writeToStdout;

    TestApplicationBase(Impl* impl) : Base(impl), out(), writeToStdout(false) {
      setOutputFolder(".");
    }

    void setOutputFolder(std::string const& folder) {
      typename Base::IO* io = Base::getIOInterface();
      io->restartReadFolder = folder;
      io->restartWriteFolder = folder;
      io->writeFolder = folder;
    }

    void print(std::string const& line) {
      if (writeToStdout) {
        std::cout << line;
        std::cout.flush();
      }
      if (out) {
        out << line;
        out.flush();
      }
    }

    void setOutputFile(std::string const& file) {
      out.close();
      out.open(file);
    }
};
