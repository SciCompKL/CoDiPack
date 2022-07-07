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
#include "applicationInterface.hpp"

/** \copydoc codi::Namespace */
namespace codi {
  namespace algorithms {

    enum class RecodingInputOutputFlags {
      InY,
      InP,
      InX,
      OutY,
      OutZ,
      MaxElement
    };
    using RecordingInputOutput = EnumBitset<RecodingInputOutputFlags>;

#define ENUM RecodingInputOutputFlags
#include "../../../misc/enumOperations.tpp"

    enum class TapeEvaluationFlags {
      G,
      F,
      P,
      MaxElement
    };
    using TapeEvaluation = EnumBitset<TapeEvaluationFlags>;

#define ENUM TapeEvaluationFlags
#include "../../../misc/enumOperations.tpp"

    enum class EvaluationInputOutputFlags {
      SetY,
      SetZ,
      GetY,
      GetP,
      GetX,
      UpdateY,
      UpdateP,
      UpdateX,
      MaxElement
    };
    using EvaluationInputOutput = EnumBitset<EvaluationInputOutputFlags>;

#define ENUM EvaluationInputOutputFlags
#include "../../../misc/enumOperations.tpp"

    template<typename T_App>
    struct AlgorithmData {
      public:
        using App = CODI_DD(T_App, CODI_T(ApplicationInterface<CODI_ANY>));
        using Type = typename App::Type;

        using Real = typename Type::Real;
        using Identifier = typename Type::Identifier;

        using RealVector = std::vector<Real>;
        using IdVector = std::vector<Identifier>;

        AlgorithmData() = default;
        AlgorithmData(App& app) : AlgorithmData() {
          this->init(app);
        }

        IdVector idInY;
        IdVector idInP;
        IdVector idInX;
        IdVector idOutY;
        IdVector idOutZ;

        RealVector realCurY;
        RealVector realNextY;

        RealVector realP;
        RealVector realX;

        void init(App& app) {
          idInY.resize(app.getSizeY());
          idInP.resize(app.getSizeP());
          idInX.resize(app.getSizeX());
          idOutY.resize(app.getSizeY());
          idOutZ.resize(app.getSizeZ());

          realCurY.resize(app.getSizeY());
          realNextY.resize(app.getSizeY());

          realP.resize(app.getSizeP());
          realX.resize(app.getSizeX());
        }
    };

    template<typename T_App>
    struct AlgorithmInterface {
      public:

        using App = CODI_DD(T_App, CODI_T(ApplicationInterface<CODI_ANY>));
        using Type = typename App::Type;

        using Real = typename Type::Real;
        using Identifier = typename Type::Identifier;
        using Tape = typename Type::Tape;

        using Data = AlgorithmData<App>;
        using Res = Residuum<Real>;

        using RealVector = typename Data::RealVector;
        using IdVector = typename Data::IdVector;

        void run(App& app);

      protected:

        void iterateUntil(App& app, int iteration) {
          while (app.getIteration() < iteration) {
            app.evaluateG();
          }
        }

        void recordTape(App& app, Data& data, TapeEvaluation evalOpt, RecordingInputOutput recOpt) {
          Tape& tape = Type::getTape();
          tape.reset();
          tape.setActive();

          if (RecodingInputOutputFlags::InY & recOpt) {
            app.iterateY(RegisterInput(data.idInY));
          } else {
            app.iterateY(clearInput);
          }
          if (RecodingInputOutputFlags::InP & recOpt) {
            app.iterateP(RegisterInput(data.idInP));
          } else {
            app.iterateP(clearInput);
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

          if (RecodingInputOutputFlags::OutZ & recOpt) {
            app.iterateZ(RegisterOutput(data.idOutZ));
          }

          tape.setPassive();
        }

        void evaluateTape(Data& data, EvaluationInputOutput operations) {
          Tape& tape = Type::getTape();

          if (EvaluationInputOutputFlags::SetY & operations) {
            setGradient(tape, data.idOutY, data.realCurY);
          }

          if (EvaluationInputOutputFlags::SetZ & operations) {
            setGradient(tape, data.idOutZ, 1.0);
          }

          tape.evaluate();

          if (EvaluationInputOutputFlags::GetY & operations) {
            getGradientAndReset(tape, data.idInY, data.realNextY);
          } else if (EvaluationInputOutputFlags::UpdateY & operations) {
            updateGradientAndReset(tape, data.idInY, data.realNextY);
          }

          if (EvaluationInputOutputFlags::GetP & operations) {
            getGradientAndReset(tape, data.idInP, data.realP);
          } else if (EvaluationInputOutputFlags::UpdateP & operations) {
            updateGradientAndReset(tape, data.idInP, data.realP);
          }

          if (EvaluationInputOutputFlags::GetX & operations) {
            getGradientAndReset(tape, data.idInX, data.realX);
          } else if (EvaluationInputOutputFlags::UpdateX & operations) {
            updateGradientAndReset(tape, data.idInX, data.realX);
          }
        }

        static void clearInput(Type& value, size_t pos) {
          Type::getTape().deactivateValue(value);
        }

        struct GetPrimal {
          public:
            RealVector& vec;
            GetPrimal(RealVector& vec) : vec(vec) {}

            void operator()(Type& value, size_t pos) {
              vec[pos] = value.getValue();
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

        static void setGradient(Tape& tape, IdVector& ids, RealVector& seed) {
          for (size_t pos = 0; pos < ids.size(); pos += 1) {
            tape.setGradient(ids[pos], seed[pos]);
          }
        }

        static void setGradient(Tape& tape, IdVector& ids, Real const& seed) {
          for (size_t pos = 0; pos < ids.size(); pos += 1) {
            tape.setGradient(ids[pos], seed);
          }
        }

        static void getGradientAndReset(Tape& tape, IdVector& ids, RealVector& value) {
          for (size_t pos = 0; pos < ids.size(); pos += 1) {
            Real& gradient = tape.gradient(ids[pos]);
            value[pos] = gradient;
            gradient = 0;
          }
        }

        static void updateGradientAndReset(Tape& tape, IdVector& ids, RealVector& value) {
          for (size_t pos = 0; pos < ids.size(); pos += 1) {
            Real& gradient = tape.gradient(ids[pos]);
            value[pos] += gradient;
            gradient = 0;
          }
        }

        static void copyFromTo(RealVector const& from, RealVector& to) {
          for (size_t pos = 0; pos < from.size(); pos += 1) {
            to[pos] = from[pos];
          }
        }
    };
  }
}
