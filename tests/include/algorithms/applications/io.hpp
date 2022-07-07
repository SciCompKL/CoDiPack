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

#include "../../utils/fileSystem.hpp"
#include "applicationBaseSettings.hpp"

template<typename T_Type, typename T_Application>
struct TestIO : public codi::algorithms::IOInterface<T_Type> {

    using Type = CODI_DD(T_Type, CODI_T(codi::LhsExpressionInterface<double, double, CODI_ANY, CODI_ANY>));
    using Application = CODI_DD(T_Application, CODI_T(codi::algorithms::ApplicationInterface<CODI_ANY>));

    using Real = typename Type::Real;

    std::map<std::string, std::vector<Real>> restartVector;
    std::map<std::string, std::pair<char*, size_t>> restartData;

    Application* app;
    TestApplicationBaseSettings& generalSettings;

    TestIO(Application* app) :
      restartVector(),
      restartData(),
      app(app),
      generalSettings(app->generalSettings)
    {}

    void writeRestartY(std::string const& fileName, std::vector<Real> const& v) {
      restartVector[fileName] = v;
    }

    void writeRestartX(std::string const& fileName, std::vector<Real> const& v) {
      restartVector[fileName] = v;
    }
    void writeRestartP(std::string const& fileName, std::vector<Real> const& v) {
      restartVector[fileName] = v;
    }

    void writeRestartData(std::string const& fileName, char* data, size_t length) {
      restartData[fileName] = std::make_pair(data, length);
    }

    void readRestartY(std::string const& fileName, std::vector<Real>& v) {
      validate(restartVector, fileName);

      v = restartVector[fileName];
    }

    void readRestartX(std::string const& fileName, std::vector<Real>& v) {
      validate(restartVector, fileName);

      v = restartVector[fileName];
    }

    void readRestartP(std::string const& fileName, std::vector<Real>& v) {
      validate(restartVector, fileName);

      v = restartVector[fileName];
    }

    void readRestartData(std::string const& fileName, char*& data, size_t& length) {
      validate(restartData, fileName);

      std::pair<char*, size_t> d = restartData[fileName];
      data = d.first;
      length = d.second;
    }

    void writeY(int iteration, std::vector<Real> const& v, codi::algorithms::OutputHints flags) {
      if(generalSettings.writeY) {
        writeVectorToFile(codi::StringUtil::format("y_%05d.txt", iteration), v, flags);
      }
    }

    void writeX(int iteration, std::vector<Real> const& v, codi::algorithms::OutputHints flags) {
      if(generalSettings.writeX) {
        writeVectorToFile(codi::StringUtil::format("x_%05d.txt", iteration), v, flags);
      }
    }

    void writeP(int iteration, std::vector<Real> const& v, codi::algorithms::OutputHints flags) {
      if(generalSettings.writeP) {
        writeVectorToFile(codi::StringUtil::format("p_%05d.txt", iteration), v, flags);
      }
    }

    void writeZ(int iteration, std::vector<Real> const& v, codi::algorithms::OutputHints flags) {
      if(generalSettings.writeZ) {
        writeVectorToFile(codi::StringUtil::format("z_%05d.txt", iteration), v, flags);
      }
    }

  private:
    template<typename D>
    void validate(std::map<std::string, D> const& map, std::string const& key) {
      if(map.end() != map.find(key)) {
        std::cerr << codi::StringUtil::format("Key '%s' not found in restart data.", key.c_str()) << std::endl;
      }
    }    

    void writeVectorToFile(std::string const& fileName, std::vector<Real> const& v, codi::algorithms::OutputHints flags) {
      if(generalSettings.onlyWriteFinal && !(codi::algorithms::OutputFlags::Final & flags)) {
        return;
      }

      FileSystem::makePath(generalSettings.outputDir.c_str());

      std::string prefix = "";
      if(codi::algorithms::OutputFlags::Primal & flags) {
        prefix = "primal_";
      } else if(codi::algorithms::OutputFlags::Derivative & flags) {
        prefix = "deriv_";
      }

      std::ofstream out(codi::StringUtil::format("%s/%s%s", generalSettings.outputDir.c_str(), prefix.c_str(), fileName.c_str()));
      for(size_t i = 0; i < v.size(); i += 1) {
        out << codi::StringUtil::format("%0.12e\n", v[i]);
      }
      out.close();
    }
};
