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

          bool pStateAvailable = app.getHints() & ApplicationHintsFlags::PStateIsAvailable;
          bool fCompAvailable = app.getHints() & ApplicationHintsFlags::FComputationIsAvailable;

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

          TapeRecordingInputOutput tapeStatus;
          TapeEvaluationInputOutput evalFlags;
          data.init(app);
          std::vector<Res> initalResY(app.getNumberOfFunctionals());

          if(fCompAvailable) {
            initializeAndComputeF(app, data, tapeStatus);
          }

          app.print(Base::formatAdjointHeader(initalResY));

          int curAdjIteration = 0;
          while (!(isFinished || isStop || isConverged)) {
            TapeRecordingInputOutput G_FLAGS = TapeRecodingInputOutputFlags::InY | TapeRecodingInputOutputFlags::OutY;
            if(!fCompAvailable) {
              G_FLAGS |= TapeRecodingInputOutputFlags::OutZ;
            }
            if (G_FLAGS != tapeStatus) {
              cpm->load(cp);
              tapeStatus = G_FLAGS;
              Base::recordTape(app, data, TapeEvaluationFlags::G, tapeStatus);
            }

            if(fCompAvailable) {
              Base::copyFromTo(yRealF, data.realNextY);
              evalFlags = TapeEvaluationInputOutputFlags::UpdateY | TapeEvaluationInputOutputFlags::SetY;
            } else {
              evalFlags = TapeEvaluationInputOutputFlags::GetY | TapeEvaluationInputOutputFlags::SetY | TapeEvaluationInputOutputFlags::SetZ;
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
                       FileOutputHintsFlags::Derivative | FileOutputHintsFlags::G | FileOutputHintsFlags::Intermediate);

            isFinished = curAdjIteration >= settings.maxIterations;
            if (1 == curAdjIteration) {
              initalResY = resY;
            } else {
              isConverged = checkConvergence(initalResY, resY);
            }
            isStop = app.isStop();
          }

          cpm->load(cp);
          tapeStatus = TapeRecodingInputOutputFlags::InX | TapeRecodingInputOutputFlags::InP | TapeRecodingInputOutputFlags::OutY;
          if(!fCompAvailable) {
            tapeStatus |= TapeRecodingInputOutputFlags::OutZ;
          }
          TapeEvaluation tapeEvaluation = TapeEvaluationFlags::G;
          if(app.getHints() & ApplicationHintsFlags::PComputationIsAvailable && !pStateAvailable) {
            tapeEvaluation |= TapeEvaluationFlags::P; // P but not iterable, then P and G need to be called at the same time.
          }
          Base::recordTape(app, data, tapeEvaluation, tapeStatus);

          if(fCompAvailable) {
            Base::copyFromTo(pRealF, data.realP);
            Base::copyFromTo(xRealF, data.realX);
            evalFlags = TapeEvaluationInputOutputFlags::SetY | TapeEvaluationInputOutputFlags::UpdateX | TapeEvaluationInputOutputFlags::UpdateP;
          } else {
            evalFlags = TapeEvaluationInputOutputFlags::SetY | TapeEvaluationInputOutputFlags::SetZ | TapeEvaluationInputOutputFlags::GetX | TapeEvaluationInputOutputFlags::GetP;
          }
          Base::evaluateTape(app, data, evalFlags);

          if(pStateAvailable) {
            Base::reverseP(app, data, TapeEvaluationInputOutputFlags::UpdateX);
          }

          io->writeY(curAdjIteration, data.realCurY, FileOutputHintsFlags::Derivative | FileOutputHintsFlags::G | FileOutputHintsFlags::Final);
          io->writeX(curAdjIteration, data.realX, FileOutputHintsFlags::Derivative | FileOutputHintsFlags::P | FileOutputHintsFlags::Final);
          if(pStateAvailable) {
            io->writeP(curAdjIteration, data.realP, FileOutputHintsFlags::Derivative | FileOutputHintsFlags::G | FileOutputHintsFlags::Final);
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

        void initializeAndComputeF(App& app, Data& data, TapeRecordingInputOutput& tapeStatus) {
          ApplicationIOInterface<Type>* io = app.getIOInterface();

          bool pStateAvailable = app.getHints() & ApplicationHintsFlags::PStateIsAvailable;

          tapeStatus = TapeRecodingInputOutputFlags::InP | TapeRecodingInputOutputFlags::InX | TapeRecodingInputOutputFlags::InY |
                       TapeRecodingInputOutputFlags::OutZ;
          Base::recordTape(app, data, TapeEvaluationFlags::F, tapeStatus);
          Base::evaluateTape(app, data, TapeEvaluationInputOutputFlags::GetP | TapeEvaluationInputOutputFlags::GetX |
                                       TapeEvaluationInputOutputFlags::GetY | TapeEvaluationInputOutputFlags::SetZ);

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

          io->writeY(0, yRealF, FileOutputHintsFlags::Derivative | FileOutputHintsFlags::F | FileOutputHintsFlags::Intermediate);
          io->writeX(0, xRealF, FileOutputHintsFlags::Derivative | FileOutputHintsFlags::F | FileOutputHintsFlags::Intermediate);
          if(pStateAvailable) {
            io->writeP(0, pRealF, FileOutputHintsFlags::Derivative | FileOutputHintsFlags::F | FileOutputHintsFlags::Intermediate);
          }
        }

        void freeF(App& app) {
          bool pStateAvailable = app.getHints() & ApplicationHintsFlags::PStateIsAvailable;

          yRealF.resize(0);
          xRealF.resize(0);
          if(pStateAvailable) {
            pRealF.resize(0);
          }
        }
    };
  }
}
