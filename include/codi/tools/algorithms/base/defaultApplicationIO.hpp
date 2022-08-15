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

#include <string>
#include <vector>

#include "../../../config.h"
#include "../../../expressions/lhsExpressionInterface.hpp"
#include "../../../misc/macros.hpp"
#include "../../../misc/stringUtil.hpp"
#include "../interfaces/applicationInterface.hpp"
#include "../interfaces/applicationIOInterface.hpp"
#include "../interfaces/fileIOInterface.hpp"

/** \copydoc codi::Namespace */
namespace codi {
  namespace algorithms {

    template<typename T_Type, typename T_WriteIO, typename T_RestartIO, typename T_Application>
    struct DefaultApplicationIO : public codi::algorithms::ApplicationIOInterface<T_Type> {
        using Type = CODI_DD(T_Type, CODI_T(LhsExpressionInterface<double, double, CODI_ANY, CODI_ANY>));
        using WriteIO = CODI_DD(T_WriteIO, FileIOInterface);
        using RestartIO = CODI_DD(T_RestartIO, FileIOInterface);
        using Application = CODI_DD(T_Application, CODI_T(codi::algorithms::ApplicationInterface<CODI_ANY>));

        using Real = RealTraits::Real<Type>;

        std::string restartWriteFolder;
        std::string restartReadFolder;
        std::string writeFolder;

        Application* app;
        WriteIO* writeIO;
        RestartIO* restartIO;

        bool outputY;
        bool outputX;
        bool outputP;
        bool outputZ;
        bool onlyWriteFinal;

        DefaultApplicationIO(Application* app, WriteIO* writeIO, RestartIO* restartIO)
            : restartWriteFolder("."),
              restartReadFolder("."),
              writeFolder("."),
              app(app),
              writeIO(writeIO),
              restartIO(restartIO),
              outputY(true),
              outputX(true),
              outputP(true),
              outputZ(true),
              onlyWriteFinal(true) {}

        void writeRestartY(std::string const& fileName, std::vector<Real> const& v) {
          writeVector(createRestartName(restartWriteFolder, fileName), restartIO, v.data(), v.size());
        }

        void writeRestartX(std::string const& fileName, std::vector<Real> const& v) {
          writeVector(createRestartName(restartWriteFolder, fileName), restartIO, v.data(), v.size());
        }
        void writeRestartP(std::string const& fileName, std::vector<Real> const& v) {
          writeVector(createRestartName(restartWriteFolder, fileName), restartIO, v.data(), v.size());
        }

        void writeRestartData(std::string const& fileName, char* data, size_t length) {
          writeVector(createRestartName(restartWriteFolder, fileName), restartIO, data, length);
        }

        void readRestartY(std::string const& fileName, std::vector<Real>& v) {
          readVector(createRestartName(restartReadFolder, fileName), restartIO, v.data(), v.size());
        }

        void readRestartX(std::string const& fileName, std::vector<Real>& v) {
          readVector(createRestartName(restartReadFolder, fileName), restartIO, v.data(), v.size());
        }

        void readRestartP(std::string const& fileName, std::vector<Real>& v) {
          readVector(createRestartName(restartReadFolder, fileName), restartIO, v.data(), v.size());
        }

        void readRestartData(std::string const& fileName, char*& data, size_t& length) {
          readVector(createRestartName(restartReadFolder, fileName), restartIO, data, length);
        }

        void writeY(int iteration, std::vector<Real> const& v, codi::algorithms::OutputHints flags, int vec) {
          if (outputY && checkFinal(flags)) {
            writeVector(createWriteName(writeFolder, "y", iteration, flags, writeIO, vec), writeIO, v.data(), v.size());
          }
        }

        void writeX(int iteration, std::vector<Real> const& v, codi::algorithms::OutputHints flags, int vec) {
          if (outputX && checkFinal(flags)) {
            writeVector(createWriteName(writeFolder, "x", iteration, flags, writeIO, vec), writeIO, v.data(), v.size());
          }
        }

        void writeP(int iteration, std::vector<Real> const& v, codi::algorithms::OutputHints flags, int vec) {
          if (outputP && checkFinal(flags)) {
            writeVector(createWriteName(writeFolder, "p", iteration, flags, writeIO, vec), writeIO, v.data(), v.size());
          }
        }

        void writeZ(int iteration, std::vector<Real> const& v, codi::algorithms::OutputHints flags, int vec) {
          if (outputZ && checkFinal(flags)) {
            writeVector(createWriteName(writeFolder, "z", iteration, flags, writeIO, vec), writeIO, v.data(), v.size());
          }
        }

      protected:

        std::string createRestartName(std::string const& folder, std::string const& name) {
          return folder + "/" + name;
        }

        template<typename IO>
        std::string createWriteName(std::string const& folder, std::string const& name, int iteration,
                                    OutputHints flags, CODI_DD(IO, FileIOInterface) * io, int vec) {
          std::string prefix = "";
          if (codi::algorithms::OutputFlags::Primal & flags) {
            prefix = "primal_";
          } else if (codi::algorithms::OutputFlags::Derivative & flags) {
            prefix = "deriv_";
          }
          if (codi::algorithms::OutputFlags::V1 & flags) {
            prefix += "v1_";
          } else if (codi::algorithms::OutputFlags::V2 & flags) {
            prefix += "v2_";
          }

          std::string suffix = "";
          if((OutputFlags::Derivative & flags) && (app->getNumberOfFunctionals() != 1)) {
            suffix += StringUtil::format("_%02d", vec);
          }

          return StringUtil::format("%s/%s%s_%05d%s.%s", folder.c_str(), prefix.c_str(), name.c_str(), iteration,
                                    suffix.c_str(), io->getFileEnding().c_str());
        }

        template<typename IO, typename T>
        void writeVector(std::string const& filename, CODI_DD(IO, FileIOInterface) * io, T const* data, size_t size) {
          size_t totalSize = sizeof(T) * size;
          typename IO::WriteHandle h = io->openWrite(filename, totalSize);
          io->write(h, data, size);
          io->closeWrite(h);
        }

        template<typename IO, typename T>
        void readVector(std::string const& filename, CODI_DD(IO, FileIOInterface) * io, T* data, size_t size) {
          typename IO::WriteHandle h = io->openRead(filename);
          io->read(h, data, size);
          io->closeRead(h);
        }

        void checkVector(std::vector<Real>& v, size_t targetSize) {
          if (v.size() != targetSize) {
            v.resize(targetSize);
          }
        }

        bool checkFinal(OutputHints flags) {
          if (onlyWriteFinal) {
            return flags & OutputFlags::Final;
          } else {
            return true;
          }
        }
    };
  }
}
