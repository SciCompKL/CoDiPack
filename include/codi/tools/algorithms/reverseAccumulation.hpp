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
#include "base/algorithmBase.hpp"

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
    struct ReverseAccumulation : public AlgorithmBase<T_App> {
      public:

        using App = CODI_DD(T_App, CODI_T(ApplicationInterface<CODI_ANY>));

        using Type = typename App::Type;

        using Base = AlgorithmBase<App>;
        using Data = typename Base::Data;
        using Res = typename Base::Res;
        using RealVector = typename Base::RealVector;

        ReverseAccumulationSettings settings;

        std::vector<RealVector> yRealF;
        std::vector<RealVector> pRealF;
        std::vector<RealVector> xRealF;

        ReverseAccumulation(ReverseAccumulationSettings settings) : settings(settings), yRealF(), pRealF(), xRealF() {}

        AlgorithmBaseSettings const* getSettings() const {
          return &settings;
        }

        void run(App& app) {
          // TODO: Handle no F available
          CheckpointManagerInterface* cpm = app.getCheckpointInterface();
          ApplicationIOInterface<Type>* io = app.getIOInterface();

          bool pStateAvailable = app.getHints() & ApplicationFlags::PStateIsAvailable;
          bool fCompAvailable = app.getHints() & ApplicationFlags::FComputationIsAvailable;

          if(!fCompAvailable) {
            app.print("Warning: Reverse accumulation without an f computation is mathematical not correct.\n");
          }

          Base::initVectorMode(app);

          Data data;
          Base::initializeApp(app, data);

          if(-1 != settings.start) {
            Base::loadClosestCheckPoint(app, settings.start);
          }

          if (-1 != settings.start && settings.start > app.getIteration()) {
            // Not yet at initial iteration
            Base::iterateUntil(app, settings.start);
          }

          CheckpointHandle* cp = cpm->create();

          bool isConverged = false;
          bool isStop = false;
          bool isFinished = false;

          RecordingInputOutput tapeStatus;
          EvaluationInputOutput evalFlags;
          data.init(app);
          std::vector<Res> initalResY(app.getNumberOfFunctionals());

          if(fCompAvailable) {
            initializeAndComputeF(app, data, tapeStatus);
          }

          app.print(Base::formatAdjointHeader(initalResY));

          int curAdjIteration = 0;
          while (!(isFinished || isStop || isConverged)) {
            RecordingInputOutput G_FLAGS = RecodingInputOutputFlags::InY | RecodingInputOutputFlags::OutY;
            if(!fCompAvailable) {
              G_FLAGS |= RecodingInputOutputFlags::OutZ;
            }
            if (G_FLAGS != tapeStatus) {
              cpm->load(cp);
              tapeStatus = G_FLAGS;
              Base::recordTape(app, data, TapeEvaluationFlags::G, tapeStatus);
            }

            if(fCompAvailable) {
              Base::copyFromTo(yRealF, data.realNextY);
              evalFlags = EvaluationInputOutputFlags::UpdateY | EvaluationInputOutputFlags::SetY;
            } else {
              evalFlags = EvaluationInputOutputFlags::GetY | EvaluationInputOutputFlags::SetY | EvaluationInputOutputFlags::SetZ;
            }
            Base::evaluateTape(app, data, evalFlags);

            std::vector<Res> resY(app.getNumberOfFunctionals());
            for(int i = 0; i < app.getNumberOfFunctionals(); i += 1) {
              resY[i] = app.residuumY(data.realCurY[i], data.realNextY[i]);
            }

            app.print(Base::formatAdjointEntry(curAdjIteration, resY));

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
          if(!fCompAvailable) {
            tapeStatus |= RecodingInputOutputFlags::OutZ;
          }
          TapeEvaluation tapeEvaluation = TapeEvaluationFlags::G;
          if(app.getHints() & ApplicationFlags::PComputationIsAvailable && !pStateAvailable) {
            tapeEvaluation |= TapeEvaluationFlags::P; // P but not iterable, then P and G need to be called at the same time.
          }
          Base::recordTape(app, data, tapeEvaluation, tapeStatus);

          if(fCompAvailable) {
            Base::copyFromTo(pRealF, data.realP);
            Base::copyFromTo(xRealF, data.realX);
            evalFlags = EvaluationInputOutputFlags::SetY | EvaluationInputOutputFlags::UpdateX | EvaluationInputOutputFlags::UpdateP;
          } else {
            evalFlags = EvaluationInputOutputFlags::SetY | EvaluationInputOutputFlags::SetZ | EvaluationInputOutputFlags::GetX | EvaluationInputOutputFlags::GetP;
          }
          Base::evaluateTape(app, data, evalFlags);

          if(pStateAvailable) {
            Base::reverseP(app, data, EvaluationInputOutputFlags::UpdateX);
          }

          io->writeY(curAdjIteration, data.realCurY, OutputFlags::Derivative | OutputFlags::G | OutputFlags::Final);
          io->writeX(curAdjIteration, data.realX, OutputFlags::Derivative | OutputFlags::P | OutputFlags::Final);
          if(pStateAvailable) {
            io->writeP(curAdjIteration, data.realP, OutputFlags::Derivative | OutputFlags::G | OutputFlags::Final);
          }

          if(fCompAvailable) {
            freeF(app);
          }

          cpm->remove(cp);
          cpm->free(cp);
        }

        bool checkConvergence(std::vector<Res> const& initial, std::vector<Res> const& cur) {
          bool allConverged = true;
          for(size_t i = 0; i < cur.size(); i += 1) {
            bool converged = false;
            if (settings.checkAbsConvergence) {
              converged = cur[i].l2 < settings.absThreshold;
            }
            if (settings.checkRelConvergence) {
              converged = cur[i].l2 < settings.relThreshold * initial[i].l2;
            }

            allConverged &= converged;
          }

          return allConverged;
        }

      private:

        void initializeAndComputeF(App& app, Data& data, RecordingInputOutput& tapeStatus) {
          ApplicationIOInterface<Type>* io = app.getIOInterface();

          bool pStateAvailable = app.getHints() & ApplicationFlags::PStateIsAvailable;

          tapeStatus = RecodingInputOutputFlags::InP | RecodingInputOutputFlags::InX | RecodingInputOutputFlags::InY |
                       RecodingInputOutputFlags::OutZ;
          Base::recordTape(app, data, TapeEvaluationFlags::F, tapeStatus);
          Base::evaluateTape(app, data, EvaluationInputOutputFlags::GetP | EvaluationInputOutputFlags::GetX |
                                       EvaluationInputOutputFlags::GetY | EvaluationInputOutputFlags::SetZ);

          yRealF.resize(app.getNumberOfFunctionals(), RealVector(app.getSizeY()));
          xRealF.resize(app.getNumberOfFunctionals(), RealVector(app.getSizeX()));
          if(pStateAvailable) {
            pRealF.resize(app.getNumberOfFunctionals(), RealVector(app.getSizeP()));
          }
          std::swap(data.realNextY, yRealF);
          std::swap(data.realX, xRealF);
          if(pStateAvailable) {
            std::swap(data.realP, pRealF);
          }

          Base::copyFromTo(yRealF, data.realCurY);  // Do first step.

          io->writeY(0, yRealF, OutputFlags::Derivative | OutputFlags::F | OutputFlags::Intermediate);
          io->writeX(0, xRealF, OutputFlags::Derivative | OutputFlags::F | OutputFlags::Intermediate);
          if(pStateAvailable) {
            io->writeP(0, pRealF, OutputFlags::Derivative | OutputFlags::F | OutputFlags::Intermediate);
          }
        }

        void freeF(App& app) {
          bool pStateAvailable = app.getHints() & ApplicationFlags::PStateIsAvailable;

          yRealF.resize(0);
          xRealF.resize(0);
          if(pStateAvailable) {
            pRealF.resize(0);
          }
        }
    };
  }
}
