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
#include "../../../misc/enumBitset.hpp"
#include "../../../misc/macros.hpp"
#include "../interfaces/applicationInterface.hpp"
#include "../interfaces/fileIOInterface.hpp"
#include "defaultApplicationIO.hpp"
#include "stateBasedCheckpointManager.hpp"

/** \copydoc codi::Namespace */
namespace codi {
  namespace algorithms {

    template<typename T_Type, typename T_Impl>
    struct DefaultApplication : public ApplicationInterface<T_Type> {
      public:

        using Type = CODI_DD(T_Type, CODI_T(LhsExpressionInterface<double, double, CODI_ANY, CODI_ANY>));
        using Impl = CODI_DD(T_Impl, CODI_T(ApplicationInterface<Type>));

        using Real = typename Type::Real;

        using CheckpointManager = StateBasedCheckpointManager<Type, BinaryFileIO, Impl>;
        using IO = DefaultApplicationIO<Type, TextFileIO, BinaryFileIO>;

      protected:

        int iteration;

        TextFileIO textIO;
        BinaryFileIO binaryIO;

        CheckpointManager cm;
        IO io;

      public:

        DefaultApplication(Impl* impl)
            : iteration(), textIO(), binaryIO(), cm("checkpoints", impl, &binaryIO), io(&textIO, &binaryIO) {
          io.restartReadFolder = "restart";
          io.restartWriteFolder = "restart";
          io.writeFolder = "output";
          io.outputY = true;
          io.outputX = true;
          io.outputP = false;
          io.outputZ = true;
          io.onlyWriteFinal = true;
        }

        virtual ~DefaultApplication() {}

        template<typename Func>
        void iterateY(Func&& func);

        template<typename Func>
        void iterateX(Func&& func);

        template<typename Func>
        void iterateP(Func&& func);

        template<typename Func>
        void iterateZ(Func&& func);

        size_t getSizeY();
        size_t getSizeX();
        size_t getSizeP();
        size_t getSizeZ();

        void evaluateG();
        void evaluateF();
        void evaluateP();

        Residuum<Real> residuumY(std::vector<Real> const& v1, std::vector<Real> const& v2) {
          return vectorBasedResiduum(v1, v2);
        }

        Residuum<Real> residuumX(std::vector<Real> const& v1, std::vector<Real> const& v2) {
          return vectorBasedResiduum(v1, v2);
        }

        Residuum<Real> residuumP(std::vector<Real> const& v1, std::vector<Real> const& v2) {
          return vectorBasedResiduum(v1, v2);
        }

        CheckpointManager* getCheckpointInterface() {
          return &cm;
        }

        IO* getIOInterface() {
          return &io;
        }

        void initialize() {}

        ApplicationHints getHints() {
          return ApplicationHints::NONE();
        }

        int getIteration() {
          return iteration;
        }

        void setIteration(int value) {
          iteration = value;
        }

        void print(std::string const& line) {
          std::cout << line;
          std::cout.flush();
        }

        bool isStop() {
          return false;
        }

      protected:

        Residuum<Real> vectorBasedResiduum(std::vector<Real> const& v1, std::vector<Real> const& v2) {
          Residuum<Real> res{};
          res.lMax = -1e300;

          codiAssert(v1.size() == v2.size());

          for (size_t i = 0; i < v1.size(); i += 1) {
            Real diff = abs(v1[i] - v2[i]);
            res.l1 += diff;
            res.l2 += diff * diff;
            if (res.lMax < diff) {
              res.lMax = diff;
              res.lMaxPos = i;
            }
          }

          res.l2 = sqrt(res.l2);

          return res;
        }
    };
  }
}
