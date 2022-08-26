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

#include "../../../config.h"
#include "../../../misc/enumBitset.hpp"
#include "../../../misc/macros.hpp"
#include "../../../traits/gradientTraits.hpp"
#include "../../helpers/customAdjointVectorHelper.hpp"
#include "../interfaces/algorithmInterface.hpp"

/** \copydoc codi::Namespace */
namespace codi {
  namespace algorithms {


    struct AlgorithmBaseSettings {
      public:

        std::string initializationTaperFolder;

        AlgorithmBaseSettings() :
          initializationTaperFolder("tapes")
        {}
    };

    template<typename T_App, typename = void>
    struct AlgorithmBase : public AlgorithmInterface<T_App> {
      public:

        using App = CODI_DD(T_App, CODI_T(ApplicationInterface<CODI_ANY>));
        using Type = typename App::Type;

        using Base = AlgorithmInterface<App>;

        using Real = typename Type::Real;
        using Identifier = typename Type::Identifier;
        using Tape = typename Type::Tape;

        using Data = AlgorithmData<App>;
        using Res = Residuum<Real>;

        using RealVector = typename Data::RealVector;
        using IdVector = typename Data::IdVector;

        using VectorAccess = VectorAccessInterface<Real, Identifier>;
        using VectorHelper = CustomAdjointVectorInterface<Type>;

        AlgorithmBase() :
          useTapeAdjoint(true),
          vectorHelper(nullptr),
          d_local(1)
        {
          setVectorMode(-1);  // Setup the default vector mode from the tape.
        }

        /// -1 for default vector mode from the tape
        void setVectorMode(int directions) {
          if(nullptr != vectorHelper) {
            delete vectorHelper;
            vectorHelper = nullptr;
          }

          int tapeVectorMode = GradientTraits::dim<typename Tape::Gradient>();
          if(-1 == directions || directions <= tapeVectorMode) {
            // Use the tape vector mode
            useTapeAdjoint = true;
            d_local = tapeVectorMode;
          } else {
            // Create a custom vector mode
            useTapeAdjoint = false;
            vectorHelper = createClosestVectorHelper(directions);
            d_local = vectorHelper->getVectorInterface()->getVectorSize();
          }
        }

        virtual AlgorithmBaseSettings const* getSettings() const = 0;

      protected:

        bool useTapeAdjoint;
        VectorHelper* vectorHelper;
        int d_local;

        using Base::iterateUntil;
        using Base::GetPrimal;

        void initVectorMode(App& app) {
          d_local = GradientTraits::dim<typename Tape::Gradient>();
        }

        void initializeApp(App& app, Data& data) {
          bool initialize = app.getHints() & ApplicationFlags::InitializationRequired;
          bool record = app.getHints() & ApplicationFlags::InitializationComputesP;
          bool pIsComputable = app.getHints() & ApplicationFlags::PComputationIsAvailable;
          bool pIsIterable = app.getHints() & ApplicationFlags::PStateIsAvailable;

          if(pIsComputable && record) {
            CODI_EXCEPTION("P can either be defined through the initialization or through the recomputation, but not "
                           "both. Manually remove either InitializationComputesP or PComputationIsAvailable from the "
                           "application hints.");
          }

          if(!pIsIterable && record) {
            CODI_EXCEPTION("P needs to be iterable if the initialization computes P.");
          }

          if(initialize || record) {

            Tape& tape = Type::getTape();

            if(record) {
              tape.reset();
              tape.setActive();
            }

            app.initialize();

            if(record) {
              data.initInitializationRecording(app);

              app.iterateX(GetId(data.idInitX));
              app.iterateP(RegisterOutput(data.idInitP));

              tape.setPassive();
              tape.swap(*data.initTape);

              if(app.getHints() & ApplicationFlags::InitializationWriteTapeToDisk) {
                data.initTape->writeToFile(getSettings()->initializationTaperFolder);
                data.initTape->deleteData();
              }
            }
          }
        }

        void loadClosestCheckPoint(App& app, int iteration) {
          CheckpointManagerInterface* cm = app.getCheckpointInterface();

          std::vector<CheckpointHandle*> checkpoints = cm->list();

          while(!checkpoints.empty() && checkpoints.back()->getIteration() > iteration) {
            cm->remove(checkpoints.back());
            checkpoints.pop_back();
          }

          if(!checkpoints.empty()) {
            cm->load(checkpoints.back());
          }

          while(!checkpoints.empty()) {
            cm->remove(checkpoints.back());
            checkpoints.pop_back();
          }
        }

        void reverseP(App& app, Data& data, EvaluationInputOutput evalXFlag) {
          if(app.getHints() & ApplicationFlags::InitializationComputesP) {

            if(app.getHints() & ApplicationFlags::InitializationWriteTapeToDisk) {
              data.initTape->readFromFile(getSettings()->initializationTaperFolder);
            }

            int d = app.getNumberOfFunctionals();

            if(!useTapeAdjoint) {
              vectorHelper->setTape(*data.initTape);
            }
            VectorAccess* access = createVectorAccess(*data.initTape);

            for(int vecPos = 0; vecPos < d; vecPos += d_local) {
              int steps = min(d - vecPos, d_local);

              setGradient(access, data.idInitP, data.realP, vecPos, steps);

              if(useTapeAdjoint) {
                data.initTape->evaluate();
              } else {
                vectorHelper->evaluate();
              }

              if (EvaluationInputOutputFlags::GetX & evalXFlag) {
                getGradientAndReset(access, data.idInitX, data.realX, vecPos, steps);
              } else if (EvaluationInputOutputFlags::UpdateX & evalXFlag) {
                updateGradientAndReset(access, data.idInitX, data.realX, vecPos, steps);
              }
            }

            deleteVectorAccess(*data.initTape, access);
            if(!useTapeAdjoint) {
              vectorHelper->setTape(Type::getTape());
            }

            if(app.getHints() & ApplicationFlags::InitializationWriteTapeToDisk) {
              data.initTape->deleteData();
            }

          } else if(app.getHints() & ApplicationFlags::PStateIsAvailable) {
            // Regular recording and reversal
            recordTape(app, data, TapeEvaluationFlags::P, RecodingInputOutputFlags::InX | RecodingInputOutputFlags::OutP);

            evaluateTape(app, data, EvaluationInputOutputFlags::SetP | evalXFlag);
          }
        }

        void recordTape(App& app, Data& data, TapeEvaluation evalOpt, RecordingInputOutput recOpt) {

          if(!useTapeAdjoint) {
            vectorHelper->deleteAdjointVector();
          }

          Tape& tape = Type::getTape();
          tape.reset();
          tape.setActive();

          if (RecodingInputOutputFlags::InY & recOpt) {
            app.iterateY(RegisterInput(data.idInY));
          } else {
            app.iterateY(clearInput);
          }
          if(app.getHints() & ApplicationFlags::PStateIsAvailable) {
            if (RecodingInputOutputFlags::InP & recOpt) {
              app.iterateP(RegisterInput(data.idInP));
            } else {
              app.iterateP(clearInput);
            }
          } else if(app.getHints() & ApplicationFlags::PComputationIsAvailable){
            evalOpt |= TapeEvaluationFlags::P; // Force the evaluation of P for clearing.
            // TODO: Force based on last tape recording.
          }
          if (RecodingInputOutputFlags::InX & recOpt) {
            app.iterateX(RegisterInput(data.idInX));
          } else {
            app.iterateX(clearInput);
          }

          if (TapeEvaluationFlags::P & evalOpt) {
            app.evaluateP();
          }
          if (TapeEvaluationFlags::G & evalOpt) {
            app.evaluateG();
          }
          if (TapeEvaluationFlags::F & evalOpt) {
            app.evaluateF();
          }

          if (RecodingInputOutputFlags::OutY & recOpt) {
            app.iterateY(RegisterOutput(data.idOutY));
          }

          if(app.getHints() & ApplicationFlags::PStateIsAvailable) {
            if (RecodingInputOutputFlags::OutP & recOpt) {
              app.iterateP(RegisterOutput(data.idOutP));
            }
          }

          if (RecodingInputOutputFlags::OutZ & recOpt) {
            app.iterateZ(RegisterOutput(data.idOutZ));
          }

          tape.setPassive();

          if(!useTapeAdjoint) {
            tape.deleteAdjointVector(); // Free memory that the tape has allocated for the adjoints.
          }
        }

        void evaluateTape(App& app, Data& data, EvaluationInputOutput operations) {
          Tape& tape = Type::getTape();

          VectorAccess* access = createVectorAccess(tape);

          int d = app.getNumberOfFunctionals();

          for(int vecPos = 0; vecPos < d; vecPos += d_local) {
            int steps = min(d - vecPos, d_local);

            if (EvaluationInputOutputFlags::SetY & operations) {
              setGradient(access, data.idOutY, data.realCurY, vecPos, steps);
            }

            if(app.getHints() & ApplicationFlags::PStateIsAvailable) {
              if (EvaluationInputOutputFlags::SetP & operations) {
                setGradient(access, data.idOutP, data.realP, vecPos, steps);
              }
            }

            if (EvaluationInputOutputFlags::SetZ & operations) {
              setGradient(access, data.idOutZ, 1.0, vecPos, steps);
            }

            if(useTapeAdjoint) {
              tape.evaluate();
            } else {
              vectorHelper->evaluate();
            }

            if (EvaluationInputOutputFlags::GetY & operations) {
              getGradientAndReset(access, data.idInY, data.realNextY, vecPos, steps);
            } else if (EvaluationInputOutputFlags::UpdateY & operations) {
              updateGradientAndReset(access, data.idInY, data.realNextY, vecPos, steps);
            }

            if(app.getHints() & ApplicationFlags::PStateIsAvailable) {
              if (EvaluationInputOutputFlags::GetP & operations) {
                getGradientAndReset(access, data.idInP, data.realP, vecPos, steps);
              } else if (EvaluationInputOutputFlags::UpdateP & operations) {
                updateGradientAndReset(access, data.idInP, data.realP, vecPos, steps);
              }
            }

            if (EvaluationInputOutputFlags::GetX & operations) {
              getGradientAndReset(access, data.idInX, data.realX, vecPos, steps);
            } else if (EvaluationInputOutputFlags::UpdateX & operations) {
              updateGradientAndReset(access, data.idInX, data.realX, vecPos, steps);
            }
          }

          deleteVectorAccess(tape, access);
        }

        static void clearInput(Type& value, size_t pos) {
          Type::getTape().deactivateValue(value);
        }

        struct GetId {
          public:
            IdVector& vec;
            GetId(IdVector& vec) : vec(vec) {}

            void operator()(Type& value, size_t pos) {
              vec[pos] = value.getIdentifier();
            }
        };

        struct RegisterInput {
          public:
            IdVector& vec;
            RegisterInput(IdVector& vec) : vec(vec) {}

            void operator()(Type& value, size_t pos) {
              Type::getTape().registerInput(value);
              vec[pos] = value.getIdentifier();
            }
        };

        struct RegisterOutput {
          public:
            IdVector& vec;
            RegisterOutput(IdVector& vec) : vec(vec) {}

            void operator()(Type& value, size_t pos) {
              Type::getTape().registerOutput(value);
              vec[pos] = value.getIdentifier();
            }
        };

        static void setGradient(VectorAccess* access, IdVector& ids, std::vector<RealVector>& seed, int vecPos, int steps) {
          RealVector vec(access->getVectorSize());
          for (size_t pos = 0; pos < ids.size(); pos += 1) {
            for (int i = 0; i < steps; i += 1) {
              vec[i] = seed[vecPos + i][pos];
            }
            access->resetAdjointVec(ids[pos]);
            access->updateAdjointVec(ids[pos], vec.data());
          }
        }

        static void setGradient(VectorAccess* access, IdVector& ids, Real const& seed, int vecPos, int steps) {
          for (int i = 0; i < steps; i += 1) {
            access->resetAdjointVec(ids[vecPos + i]);
            access->updateAdjoint(ids[vecPos + i], i, seed);
          }

        }

        static void getGradientAndReset(VectorAccess* access, IdVector& ids, std::vector<RealVector>& value, int vecPos, int steps) {
          RealVector vec(access->getVectorSize());

          for (size_t pos = 0; pos < ids.size(); pos += 1) {
            access->getAdjointVec(ids[pos], vec.data());
            access->resetAdjointVec(ids[pos]);
            for (int i = 0; i < steps; i += 1) {
              value[vecPos + i][pos] = vec[i];
            }
          }
        }

        static void updateGradientAndReset(VectorAccess* access, IdVector& ids, std::vector<RealVector>& value, int vecPos, int steps) {
          RealVector vec(access->getVectorSize());

          for (size_t pos = 0; pos < ids.size(); pos += 1) {
            access->getAdjointVec(ids[pos], vec.data());
            access->resetAdjointVec(ids[pos]);
            for (int i = 0; i < steps; i += 1) {
              value[vecPos + i][pos] += vec[i];
            }
          }
        }

        static void copyFromTo(std::vector<RealVector> const& from, std::vector<RealVector>& to) {
          for(size_t i = 0; i < to.size(); i += 1) {
            std::copy(from.begin(), from.end(), to.begin());
          }
        }


        VectorAccess* createVectorAccess(Tape& tape) {
          if(useTapeAdjoint) {
            return tape.createVectorAccess();
          } else {
            return vectorHelper->getVectorInterface();
          }
        }

        void deleteVectorAccess(Tape& tape, VectorAccess* access) {
          if(useTapeAdjoint) {
            tape.deleteVectorAccess(access);
          } else {
            // Vector helper interfaces do not need to be freed.
          }
        }

        virtual VectorHelper* createClosestVectorHelper(int directions) {
          if(directions <= 1) {
            return createVectorHelper<1>();
          } else if(directions <= 2) {
            return createVectorHelper<2>();
          } else if(directions <= 4) {
            return createVectorHelper<4>();
          } else if(directions <= 8) {
            return createVectorHelper<8>();
          } else if(directions <= 12) {
            return createVectorHelper<12>();
          } else {
            return createVectorHelper<16>();
          }
        }

        template<int dim>
        VectorHelper* createVectorHelper() {
          return new CustomAdjointVectorHelper<Type, Direction<Real, dim>>();
        }

        std::string formatAdjointHeader(std::vector<Res> res) {
          int vectorDirections = res.size();

          std::string out = "Iter";
          for(int i = 0; i < vectorDirections; i += 1) {
            std::string prefix = StringUtil::format("V%02d_Adj", i);

            if(1 == vectorDirections) {
              prefix = "Adj";
            }
            out += " " + res[i].formatHeader(prefix);
          }

          out += "\n";

          return out;
        }

        std::string formatAdjointEntry(int adjIteration, std::vector<Res> const& resY, int width = 6) {
          std::string out = StringUtil::format("%d", adjIteration);
          for(size_t i = 0; i < resY.size(); i += 1) {
            out += " " + resY[i].formatEntry(width);
          }

          out += "\n";

          return out;
        }
    };
  }
}
