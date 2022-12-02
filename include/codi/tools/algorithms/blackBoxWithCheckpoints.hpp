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

    struct BlackBoxWithCheckpointsSettings : public AlgorithmBaseSettings {
        int start;        ///< Start iteration, -1 for initial iteration.
        int end;          ///< End iteration, -1 for last checkpoint.
        int adjointStep;  ///< Increment for the adjoint step. Usually 1.

        bool checkPrimalConvergence; ///< If a primal converged state can overwrite the end setting.

        bool verbose;       ///< Write algorithm actions.

        BlackBoxWithCheckpointsSettings()
            : start(0),
              end(0),
              adjointStep(1),
              checkPrimalConvergence(false),
              verbose(false) {}
    };

    template<typename T_App>
    struct BlackBoxWithCheckpoints : public AlgorithmBase<T_App> {
      public:

        using App = CODI_DD(T_App, CODI_T(ApplicationInterface<CODI_ANY>));

        using Type = typename App::Type;

        using Base = AlgorithmBase<App>;
        using Data = typename Base::Data;
        using Res = typename Base::Res;
        using RealVector = typename Base::RealVector;

        BlackBoxWithCheckpointsSettings settings;

        BlackBoxWithCheckpoints() : settings() {}
        BlackBoxWithCheckpoints(BlackBoxWithCheckpointsSettings settings) : settings(settings) {}

        AlgorithmBaseSettings const* getSettings() const {
          return &settings;
        }

        inline int prepareCheckpointsAtEnd(CheckpointManagerInterface* cpm, std::vector<CheckpointHandle*>& checkpoints, bool fAvailable)
        {
          int curAdjIteration;
          if(-1 == settings.end) {
            curAdjIteration = checkpoints.back()->getIteration();

            if(!fAvailable) {
              curAdjIteration += settings.adjointStep; // Increase iteration by one since we do one additional G iteration to compute f, so that the last checkpoint is used.
            }

          } else {
            curAdjIteration = settings.end;
            if(!fAvailable) {
              curAdjIteration -= settings.adjointStep; // Decrease iteration by one since compute f with the last G iteration.
            }
            while(curAdjIteration < checkpoints.back()->getIteration()) {
              cpm->free(checkpoints.back());
              checkpoints.pop_back();
            }
          }

          return curAdjIteration;
        }

        inline void prepareCheckpointsAtFront(CheckpointManagerInterface* cpm, std::vector<CheckpointHandle*>& checkpoints) {
          if(0 != checkpoints.size()) {
            int removeFront = 0;
            while(settings.start > checkpoints.front()->getIteration()) {
              removeFront += 1;
            }
            if(removeFront != 0 && settings.start != checkpoints[removeFront]->getIteration()) {
              removeFront -= 1; // Start and last checkpoint are not the same skip one removal
            }
            for(int i = 0; i < removeFront; i += 1) {
              cpm->free(checkpoints[i]);
            }
            checkpoints.erase(checkpoints.begin(), checkpoints.begin() + removeFront);
          }

          // Check if we have a checkpoint at the start or before
          if(0 == checkpoints.size() || settings.start < checkpoints.front()->getIteration()) {
            CheckpointHandle* startCheck = cpm->create();
            if(settings.start < startCheck->getIteration()) {
              CODI_EXCEPTION("Can not reach start iteration with available checkpoints.");
            }

            checkpoints.insert(checkpoints.begin(), startCheck);
          }
        }

        inline void iterateUntilWithCheckpoints(App& app, int until, std::vector<CheckpointHandle*>& checkpoints,
                                                CheckpointManagerInterface* cpm, bool checkPrimalConvergence) {
          int start = app.getIteration();
          bool continueRunning = true;

          if(settings.verbose) { app.print(StringUtil::format("Iterating from %d to %d.\n", start, until)); }

          for(int curPos = start; continueRunning && curPos < until; curPos += settings.adjointStep) {
            if(start != curPos) {
              if(settings.verbose) { app.print(StringUtil::format("Creating checkpoint at %d.\n", curPos)); }

              checkpoints.push_back(cpm->create());
            }
            app.evaluateG();

            if(checkPrimalConvergence) { continueRunning &= !app.isConverged(); }
            continueRunning &= !app.isStop();
          }
        }

        inline void loadAndPopCheckpoint(CheckpointManagerInterface* cpm, std::vector<CheckpointHandle*>& checkpoints) {
          cpm->load(checkpoints.back());

          if(1 != checkpoints.size()) { // Check for last iteration.
            popCheckpoint(cpm, checkpoints);
          }
        }

        inline void popCheckpoint(CheckpointManagerInterface* cpm, std::vector<CheckpointHandle*>& checkpoints) {
          cpm->remove(checkpoints.back());
          cpm->free(checkpoints.back());
          checkpoints.pop_back();
        }

        inline void clearCheckpoints(CheckpointManagerInterface* cpm, std::vector<CheckpointHandle*>& checkpoints) {
          while(!checkpoints.empty()) {
            popCheckpoint(cpm, checkpoints);
          }
        }

        void run(App& app) {
          CheckpointManagerInterface* cpm = app.getCheckpointInterface();
          ApplicationIOInterface<Type>* io = app.getIOInterface();

          bool pStateAvailable = app.getHints() & ApplicationHintsFlags::PStateIsAvailable;
          bool fAvailable = app.getHints() & ApplicationHintsFlags::FComputationIsAvailable;

          Base::initVectorMode(app);

          Data data;
          Base::initializeApp(app, data);

          std::vector<CheckpointHandle*> checkpoints = cpm->list();

          prepareCheckpointsAtFront(cpm, checkpoints);
          int curAdjIteration = prepareCheckpointsAtEnd(cpm, checkpoints, fAvailable);

          if(settings.verbose) {
            app.print(StringUtil::format("Checkpoints avail: %d, first: %d, last: %d\n",
                checkpoints.size(),
                checkpoints.front()->getIteration(),
                checkpoints.back()->getIteration()));
          }

          bool isFirst = true;
          bool isStop = false;
          bool isFinished = false;

          TapeRecordingInputOutput tapeStatus;
          data.init(app);
          std::vector<Res> initalResY(app.getNumberOfFunctionals());
          std::vector<Res> resY(app.getNumberOfFunctionals());

          cpm->load(checkpoints.back());
          if(curAdjIteration != app.getIteration()) {
            // Iterates until app.getIteration is at curAdjIteration. No Checkpoint is written for curAdjIteration,
            // therefore we do not need to delete the checkpoint.
            iterateUntilWithCheckpoints(app, curAdjIteration, checkpoints, cpm, settings.checkPrimalConvergence);
            // Check app state
            if(app.isStop()) {
              // Stop algorithm
              clearCheckpoints(cpm, checkpoints);
              return;
            } else if(settings.checkPrimalConvergence && app.isConverged()) {
              // Readjust end to primal convergence
              settings.end = app.getIteration();
              curAdjIteration = prepareCheckpointsAtEnd(cpm, checkpoints, fAvailable);
              if(!fAvailable) {
                // In order to compute f at the current position, we have to load the last checkpoint again.
                loadAndPopCheckpoint(cpm, checkpoints);
              }
            }
          } else {
            // Remove the loaded checkpoint.
            popCheckpoint(cpm, checkpoints);
          }

          if(settings.verbose) { app.print(StringUtil::format("Computing adjoint of f at %d.\n", app.getIteration())); }

          tapeStatus = TapeRecodingInputOutputFlags::InP | TapeRecodingInputOutputFlags::InX | TapeRecodingInputOutputFlags::InY |
                       TapeRecodingInputOutputFlags::OutZ;

          if(fAvailable) {
            Base::recordTape(app, data, TapeEvaluationFlags::F, tapeStatus);
          } else {
            Base::recordTape(app, data, TapeEvaluationFlags::G, tapeStatus);
          }
          Base::evaluateTape(app, data, TapeEvaluationInputOutputFlags::GetP | TapeEvaluationInputOutputFlags::GetX |
                                       TapeEvaluationInputOutputFlags::GetY | TapeEvaluationInputOutputFlags::SetZ);

          io->writeY(curAdjIteration, data.realNextY, FileOutputHintsFlags::Derivative | FileOutputHintsFlags::F | FileOutputHintsFlags::Intermediate);
          io->writeX(curAdjIteration, data.realX, FileOutputHintsFlags::Derivative | FileOutputHintsFlags::F | FileOutputHintsFlags::Intermediate);
          if(pStateAvailable) {
            io->writeP(curAdjIteration, data.realP, FileOutputHintsFlags::Derivative | FileOutputHintsFlags::F | FileOutputHintsFlags::Intermediate);
          }
          std::swap(data.realNextY, data.realCurY);
          curAdjIteration -= settings.adjointStep;

          if(checkpoints.back()->getIteration() == curAdjIteration) {
            loadAndPopCheckpoint(cpm, checkpoints);
          } else {
            cpm->load(checkpoints.back());
          }

          if(settings.verbose) { app.print(StringUtil::format("Starting main loop.\n")); }
          app.print(Base::formatAdjointHeader(initalResY));
          while(!(isFinished || isStop)) {
            if(app.getIteration() == curAdjIteration) {
              if(settings.verbose) { app.print(StringUtil::format("Computing adjoint of G at %d.\n", app.getIteration())); }
              tapeStatus = TapeRecodingInputOutputFlags::InP | TapeRecodingInputOutputFlags::InX | TapeRecodingInputOutputFlags::InY |
                           TapeRecodingInputOutputFlags::OutY;
              TapeEvaluation tapeEval = TapeEvaluationFlags::G;
              if(!pStateAvailable) {
                tapeEval |= TapeEvaluationFlags::P;
              }
              Base::recordTape(app, data, TapeEvaluationFlags::G | TapeEvaluationFlags::P, tapeStatus);

              Base::evaluateTape(app, data, TapeEvaluationInputOutputFlags::UpdateP | TapeEvaluationInputOutputFlags::UpdateX |
                                            TapeEvaluationInputOutputFlags::GetY | TapeEvaluationInputOutputFlags::SetY);


              for(int i = 0; i < app.getNumberOfFunctionals(); i += 1) {
                resY[i] = app.residuumY(data.realCurY[i], data.realNextY[i]);
              }
              app.print(Base::formatAdjointEntry(curAdjIteration, resY));

              io->writeY(curAdjIteration, data.realNextY,
                         FileOutputHintsFlags::Derivative | FileOutputHintsFlags::G | FileOutputHintsFlags::Intermediate);
              io->writeX(curAdjIteration, data.realX,
                         FileOutputHintsFlags::Derivative | FileOutputHintsFlags::G | FileOutputHintsFlags::Intermediate);
              if(pStateAvailable) {
                io->writeP(curAdjIteration, data.realP,
                           FileOutputHintsFlags::Derivative | FileOutputHintsFlags::G | FileOutputHintsFlags::Intermediate);
              }

              // Prepare next iteration
              std::swap(data.realNextY, data.realCurY);
              isFinished = curAdjIteration == settings.start;
              curAdjIteration -= settings.adjointStep;
              if (isFirst) {
                initalResY = resY;
                isFirst = false;
              }
              isStop = app.isStop();

              loadAndPopCheckpoint(cpm, checkpoints);
            } else {
              iterateUntilWithCheckpoints(app, curAdjIteration, checkpoints, cpm, false);
            }
          }

          if(pStateAvailable) {
            if(settings.verbose) { app.print(StringUtil::format("Computing adjoint of P at %d.\n", curAdjIteration)); }
            Base::reverseP(app, data, TapeEvaluationInputOutputFlags::UpdateX);
          }

          io->writeY(settings.start, data.realCurY, FileOutputHintsFlags::Derivative | FileOutputHintsFlags::G | FileOutputHintsFlags::Final);
          io->writeX(settings.start, data.realX, FileOutputHintsFlags::Derivative | FileOutputHintsFlags::P | FileOutputHintsFlags::Final);
          if(pStateAvailable) {
            io->writeP(settings.start, data.realP, FileOutputHintsFlags::Derivative | FileOutputHintsFlags::G | FileOutputHintsFlags::Final);
          }

          if(settings.verbose) { app.print(StringUtil::format("Finished.\n")); }
        }
    };
  }
}
