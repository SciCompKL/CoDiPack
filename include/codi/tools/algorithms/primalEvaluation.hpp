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

    struct PrimalEvaluationSettings {
        int maxIterations;  ///< Maximum number of adjoint iterations.

        bool checkAbsConvergence;
        bool checkRelConvergence;

        double absThreshold;
        double relThreshold;

        bool writeCheckpoints;
        bool writeFinalCheckpoint;
        int checkpointsInterleave;

        PrimalEvaluationSettings()
            : maxIterations(1000),
              checkAbsConvergence(true),
              checkRelConvergence(true),
              absThreshold(1e-12),
              relThreshold(1e-6),
              writeCheckpoints(false),
              writeFinalCheckpoint(false),
              checkpointsInterleave(10) {}
    };

    template<typename T_App>
    struct PrimalEvaluation : public AlgorithmInterface<T_App> {
      public:

        using App = CODI_DD(T_App, CODI_T(ApplicationInterface<CODI_ANY>));

        using Type = typename App::Type;

        using Base = AlgorithmInterface<App>;
        using Res = typename App::Res;
        using RealVector = typename Base::RealVector;

        PrimalEvaluationSettings settings;

        PrimalEvaluation() : settings() {}
        PrimalEvaluation(PrimalEvaluationSettings settings) : settings(settings) {}

        PrimalEvaluationSettings const* getSettings() const {
          return &settings;
        }

        void run(App& app) {
          ApplicationIOInterface<Type>* io = app.getIOInterface();

          bool isConverged = false;
          bool isStop = false;
          bool isFinished = false;

          RealVector yCur(app.getSizeY());
          RealVector yNext(app.getSizeY());
          Res initalResY;

          app.print(formatHeader(initalResY));

          app.initialize();
          app.evaluateP();

          app.iterateY(typename Base::GetPrimal(yCur));

          while (!(isFinished || isStop || isConverged)) {

            writeCheckpoint(app);

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

            io->writeY(app.getIteration(), yCur,
                       OutputFlags::Primal | OutputFlags::G |
                           ((isFinished | isConverged | isStop) ? OutputFlags::Final : OutputFlags::Intermediate));
          }

          writeCheckpoint(app, true);

          app.evaluateF();

          RealVector z(app.getSizeZ());
          app.iterateZ(typename Base::GetPrimal(z));

          io->writeZ(app.getIteration(), z, OutputFlags::Primal | OutputFlags::F | OutputFlags::Final);
        }

        void writeCheckpoint(App& app, bool final = false) {
          bool write = false;
          if(settings.writeCheckpoints && (0 == app.getIteration() % settings.checkpointsInterleave)) {
            write = true;
          } else if(settings.writeFinalCheckpoint && final) {
            write = true;
          }

          if(write) {
            CheckpointManagerInterface* cpm = app.getCheckpointInterface();
            CheckpointHandle* cp = cpm->create();
            cpm->write(cp);
            cpm->free(cp);
          }
        }

        std::string formatHeader(Res resY) {
          return "Iter " + resY.formatHeader("") + "\n";
        }

        std::string formatEntry(int iteration, Res resY) {
          return StringUtil::format("%d ", iteration) + resY.formatEntry() + "\n";
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
