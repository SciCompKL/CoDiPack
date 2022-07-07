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

#include <vector>

#include "../../config.h"
#include "../../expressions/lhsExpressionInterface.hpp"
#include "../../misc/macros.hpp"
#include "../../misc/stringUtil.hpp"
#include "interfaces/algorithmInterface.hpp"

/** \copydoc codi::Namespace */
namespace codi {
  namespace algorithms {

    struct BlackBoxSettings {
        int maxIterations;  ///< Maximum number of adjoint iterations.

        bool checkAbsConvergence;
        bool checkRelConvergence;

        double absThreshold;
        double relThreshold;

        BlackBoxSettings()
            : maxIterations(1000),
              checkAbsConvergence(true),
              checkRelConvergence(true),
              absThreshold(1e-12),
              relThreshold(1e-6) {}
    };

    template<typename T_App>
    struct BlackBox : public AlgorithmInterface<T_App> {
      public:

        using App = CODI_DD(T_App, CODI_T(ApplicationInterface<CODI_ANY>));

        using Type = typename App::Type;
        using Tape = typename Type::Tape;

        using Base = AlgorithmInterface<App>;
        using Data = typename Base::Data;
        using Res = typename Base::Res;
        using RealVector = typename Base::RealVector;
        using IdVector = typename Base::IdVector;

        BlackBoxSettings settings;

        BlackBox(BlackBoxSettings settings) : settings(settings) {}

        void run(App& app) {
          IOInterface<Type>* io = app.getIOInterface();

          bool isConverged = false;
          bool isStop = false;
          bool isFinished = false;

          RealVector yCur(app.getSizeY());
          RealVector yNext(app.getSizeY());
          IdVector idX(app.getSizeX());
          IdVector idZ(app.getSizeZ());
          RealVector gradX(app.getSizeX());
          Res initalResY;

          app.print(formatHeader());

          app.initialize();

          Tape& tape = Type::getTape();
          tape.setActive();
          app.iterateX(typename Base::RegisterInput(idX));

          app.evaluateP();

          app.iterateY(typename Base::GetPrimal(yCur));

          while (!(isFinished || isStop || isConverged)) {
            app.evaluateG();

            app.iterateY(typename Base::GetPrimal(yNext));

            Res resY = app.residuumY(yCur, yNext);

            app.print(formatEntry(app.getIteration(), resY));

            // Prepare next iteration
            std::swap(yCur, yNext);

            isFinished = app.getIteration() >= settings.maxIterations;
            if (1 == app.getIteration()) {
              initalResY = resY;
            } else {
              isConverged = checkConvergence(initalResY, resY);
            }
            isStop = app.isStop();
          }

          app.evaluateF();
          app.iterateZ(typename Base::RegisterOutput(idZ));

          tape.setPassive();

          Base::setGradient(tape, idZ, 1.0);

          tape.evaluate();

          Base::getGradientAndReset(tape, idX, gradX);

          io->writeX(0, gradX, OutputFlags::Final | OutputFlags::Derivative | OutputFlags::F);
        }

        std::string formatHeader() {
          return StringUtil::format("Iter Y_L1 Y_L2 Y_LMax Y_LMaxPos\n");
        }

        std::string formatEntry(int iteration, Res resY) {
          return StringUtil::format("%d %0.6e %0.6e %0.6e %d\n", iteration, resY.l1, resY.l2, resY.lMax, resY.lMaxPos);
        }

        bool checkConvergence(Res initial, Res cur) {
          bool converged = false;
          if (settings.checkAbsConvergence) {
            converged = cur.l2 < settings.absThreshold;
          }
          if (settings.checkRelConvergence) {
            converged = cur.l2 < settings.relThreshold * initial.l2;
          }

          return converged;
        }
    };
  }
}
