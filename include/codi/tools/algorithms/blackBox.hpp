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
#include "base/algorithmBase.hpp"
#include "interfaces/algorithmInterface.hpp"
#include "tools/reverseTapeOutput.hpp"

/** \copydoc codi::Namespace */
namespace codi {
  namespace algorithms {

    struct BlackBoxSettings : public AlgorithmBaseSettings {
        int maxIterations;  ///< Maximum number of adjoint iterations.

        bool outputPrimalConvergence;
        bool checkPrimalConvergence;

        bool intermediateReverseResultsOutput;

        BlackBoxSettings()
            : maxIterations(1000),
              outputPrimalConvergence(true),
              checkPrimalConvergence(true),
              intermediateReverseResultsOutput(false) {}
    };

    template<typename T_App>
    struct BlackBox : public AlgorithmBase<T_App> {
      public:

        using App = CODI_DD(T_App, CODI_T(ApplicationInterface<CODI_LHS_EXPRESSION_PROXY>));

        using Type = typename App::Type;
        using Tape = typename Type::Tape;

        using Base = AlgorithmBase<T_App>;
        using Data = typename Base::Data;
        using Res = typename Base::Res;
        using RealVector = typename Base::RealVector;
        using IdVector = typename Base::IdVector;

        BlackBoxSettings settings;

        BlackBox(BlackBoxSettings settings) : settings(settings) {}

        AlgorithmBaseSettings const* getSettings() const {
          return &settings;
        }

        void run(App& app) {
          ApplicationIOInterface<Type>* io = app.getIOInterface();

          bool continueRunning = true;

          app.initialize();

          Base::initVectorMode(app);

          RealVector yCur(0);
          RealVector yNext(0);
          if(settings.outputPrimalConvergence) {
            yCur.resize(app.getSizeY());
            yNext.resize(app.getSizeY());

            app.print(formatHeader());
          }
          IdVector idX(app.getSizeX());
          IdVector idZ(app.getSizeZ());
          std::vector<RealVector> gradX(Base::d_local, RealVector(app.getSizeX()));

          Tape& tape = Type::getTape();
          tape.setActive();
          app.iterateX(typename Base::RegisterInput(idX));

          app.evaluateP();

          if(settings.outputPrimalConvergence) {
            app.iterateY(typename Base::GetPrimal(yCur));
          }

          while (continueRunning) {
            app.evaluateG();

            if(settings.outputPrimalConvergence) {
              app.iterateY(typename Base::GetPrimal(yNext));

              Res resY = app.residuumY(yCur, yNext);

              app.print(formatEntry(app.getIteration(), resY));

              // Prepare next iteration
              std::swap(yCur, yNext);
            }

            if(settings.intermediateReverseResultsOutput) { addDebugOutput(app); }

            if(settings.checkPrimalConvergence) { continueRunning &= !app.isConverged(); }
            continueRunning &= app.getIteration() < settings.maxIterations;
            continueRunning &= !app.isStop();
          }

          app.evaluateF();
          app.iterateZ(typename Base::RegisterOutput(idZ));

          tape.setPassive();

          int d = app.getNumberOfFunctionals();

          typename Base::VectorAccess* access = Base::createVectorAccess(tape);

          for(int vecPos = 0; vecPos < d; vecPos += Base::d_local) {
            int steps = min(d - vecPos, Base::d_local);

            Base::setGradient(access, idZ, 1.0, vecPos, steps);

            if(Base::useTapeAdjoint) {
              tape.evaluate();
            } else {
              Base::vectorHelper->evaluate();
            }

            Base::getGradientAndReset(access, idX, gradX, vecPos, steps);

            io->writeX(0, gradX, FileOutputHintsFlags::Final | FileOutputHintsFlags::Derivative | FileOutputHintsFlags::F, vecPos);
          }

          Base::deleteVectorAccess(tape, access);
        }

        std::string formatHeader() {
          return StringUtil::format("Iter Y_L1 Y_L2 Y_LMax Y_LMaxPos\n");
        }

        std::string formatEntry(int iteration, Res resY) {
          return StringUtil::format("%d %0.6e %0.6e %0.6e %d\n", iteration, resY.l1, resY.l2, resY.lMax, resY.lMaxPos);
        }

      private:

        void addDebugOutput(App& app) {
          FileOutputHints hints = FileOutputHintsFlags::Intermediate | FileOutputHintsFlags::Derivative | FileOutputHintsFlags::G;

          IdVector idY(app.getSizeY());
          app.iterateY(typename Base::GetId(idY));
          ReverseTapeOutput<App>::addReverseOutput(app, idY, FileOutputType::Y, hints);

          IdVector idX(app.getSizeX());
          app.iterateX(typename Base::GetId(idX));
          ReverseTapeOutput<App>::addReverseOutput(app, idX, FileOutputType::X, hints);

          if(0 != app.getSizeP()) {
            IdVector idP(app.getSizeP());
            app.iterateP(typename Base::GetId(idP));
            ReverseTapeOutput<App>::addReverseOutput(app, idP, FileOutputType::P, hints);
          }
        }
    };
  }
}
