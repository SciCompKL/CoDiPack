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

    struct ReverseAccumulationSettings : public AlgorithmBaseSettings {
        int start;          ///< Start iteration, -1 for current position.
        int maxIterations;  ///< Maximum number of adjoint iterations.

        bool checkAbsConvergence;
        bool checkRelConvergence;

        double absThreshold;
        double relThreshold;

        ReverseAccumulationSettings()
            : start(-1),
              maxIterations(1000),
              checkAbsConvergence(true),
              checkRelConvergence(false),
              absThreshold(1e-12),
              relThreshold(1e-6) {}
    };

    template<typename T_App>
    struct ReverseAccumulation : public AlgorithmInterface<T_App> {
      public:

        using App = CODI_DD(T_App, CODI_T(ApplicationInterface<CODI_ANY>));

        using Type = typename App::Type;

        using Base = AlgorithmInterface<App>;
        using Data = typename Base::Data;
        using Res = typename Base::Res;
        using RealVector = typename Base::RealVector;

        ReverseAccumulationSettings settings;

        ReverseAccumulation(ReverseAccumulationSettings settings) : settings(settings) {}

        AlgorithmBaseSettings const* getSettings() const {
          return &settings;
        }

        void run(App& app) {
          CheckpointManagerInterface* cpm = app.getCheckpointInterface();
          ApplicationIOInterface<Type>* io = app.getIOInterface();

          Data data;
          Base::initializeApp(app, data);

          if (-1 != settings.start && settings.start > app.getIteration()) {
            // Not yet at initial iteration
            Base::iterateUntil(app, settings.start);
          }

          CheckpointHandle* cp = cpm->create();

          bool isConverged = false;
          bool isStop = false;
          bool isFinished = false;

          RecordingInputOutput tapeStatus;
          data.init(app);
          Res initalResY;

          tapeStatus = RecodingInputOutputFlags::InP | RecodingInputOutputFlags::InX | RecodingInputOutputFlags::InY |
                       RecodingInputOutputFlags::OutZ;
          Base::recordTape(app, data, TapeEvaluationFlags::F, tapeStatus);
          Base::evaluateTape(data, EvaluationInputOutputFlags::GetP | EvaluationInputOutputFlags::GetX |
                                       EvaluationInputOutputFlags::GetY | EvaluationInputOutputFlags::SetZ);

          RealVector yRealF(app.getSizeY());
          RealVector pRealF(app.getSizeP());
          RealVector xRealF(app.getSizeX());
          std::swap(data.realNextY, yRealF);
          std::swap(data.realP, pRealF);
          std::swap(data.realX, xRealF);

          Base::copyFromTo(yRealF, data.realCurY);  // Do first step.

          app.print(formatHeader());

          io->writeY(0, yRealF, OutputFlags::Derivative | OutputFlags::F | OutputFlags::Intermediate);
          io->writeX(0, xRealF, OutputFlags::Derivative | OutputFlags::F | OutputFlags::Intermediate);
          io->writeP(0, pRealF, OutputFlags::Derivative | OutputFlags::F | OutputFlags::Intermediate);

          int curAdjIteration = 0;
          while (!(isFinished || isStop || isConverged)) {
            RecordingInputOutput const G_FLAGS = RecodingInputOutputFlags::InY | RecodingInputOutputFlags::OutY;
            if (G_FLAGS != tapeStatus) {
              cpm->load(cp);
              tapeStatus = G_FLAGS;
              Base::recordTape(app, data, TapeEvaluationFlags::G, tapeStatus);
            }

            Base::copyFromTo(yRealF, data.realNextY);
            Base::evaluateTape(data, EvaluationInputOutputFlags::UpdateY | EvaluationInputOutputFlags::SetY);
            Type::getTape().clearAdjoints();  // TODO: Remove

            Res resY = app.residuumY(data.realCurY, data.realNextY);

            app.print(formatEntry(curAdjIteration, resY));

            // Prepare next iteration
            std::swap(data.realCurY, data.realNextY);
            curAdjIteration += 1;
            io->writeY(curAdjIteration, data.realCurY,
                       OutputFlags::Derivative | OutputFlags::G | OutputFlags::Intermediate);

            isFinished = curAdjIteration >= settings.maxIterations;
            if (1 == curAdjIteration) {
              initalResY = resY;
            } else {
              isConverged = checkConvergence(initalResY, resY);
            }
            isStop = app.isStop();
          }

          cpm->load(cp);
          tapeStatus = RecodingInputOutputFlags::InX | RecodingInputOutputFlags::InP | RecodingInputOutputFlags::OutY;
          Base::recordTape(app, data, TapeEvaluationFlags::G, tapeStatus);

          Base::copyFromTo(pRealF, data.realP);
          Base::copyFromTo(xRealF, data.realX);
          Base::evaluateTape(data, EvaluationInputOutputFlags::SetY | EvaluationInputOutputFlags::UpdateX | EvaluationInputOutputFlags::UpdateP);

          Base::reverseP(app, data, EvaluationInputOutputFlags::UpdateX);

          io->writeY(curAdjIteration, data.realCurY, OutputFlags::Derivative | OutputFlags::G | OutputFlags::Final);
          io->writeP(curAdjIteration, data.realP, OutputFlags::Derivative | OutputFlags::G | OutputFlags::Final);
          io->writeX(curAdjIteration, data.realX, OutputFlags::Derivative | OutputFlags::P | OutputFlags::Final);

          cpm->remove(cp);
        }

        std::string formatHeader() {
          return StringUtil::format("Iter AdjY_L1 AdjY_L2 AdjY_LMax AdjY_LMaxPos\n");
        }
        std::string formatEntry(int adjIteration, Res resY) {
          return StringUtil::format("%d %0.6e %0.6e %0.6e %d\n", adjIteration, resY.l1, resY.l2, resY.lMax,
                                    resY.lMaxPos);
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

        EvaluationInputOutput getVectorOperations(RecordingInputOutput tapeStatus, EvaluationInputOutput vectorStatus) {
          if (!(RecodingInputOutputFlags::InY & tapeStatus)) {
            // Blank out y operations
            vectorStatus.reset(EvaluationInputOutputFlags::GetY);
            vectorStatus.reset(EvaluationInputOutputFlags::UpdateY);
          }

          if (!(RecodingInputOutputFlags::InP & tapeStatus)) {
            // Blank out p operations
            vectorStatus.reset(EvaluationInputOutputFlags::GetP);
            vectorStatus.reset(EvaluationInputOutputFlags::UpdateP);
          }

          if (!(RecodingInputOutputFlags::InX & tapeStatus)) {
            // Blank out x operations
            vectorStatus.reset(EvaluationInputOutputFlags::GetX);
            vectorStatus.reset(EvaluationInputOutputFlags::UpdateX);
          }
          return vectorStatus;
        }
    };
  }
}
