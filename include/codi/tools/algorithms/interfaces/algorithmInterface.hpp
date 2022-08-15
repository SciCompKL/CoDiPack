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
#include "applicationInterface.hpp"

/** \copydoc codi::Namespace */
namespace codi {
  namespace algorithms {

    enum class RecodingInputOutputFlags {
      InY,
      InP,
      InX,
      OutY,
      OutP,
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
      SetP,
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
        using Tape = typename Type::Tape;

        using Real = typename Type::Real;
        using Identifier = typename Type::Identifier;

        using RealVector = std::vector<Real>;
        using IdVector = std::vector<Identifier>;

        AlgorithmData() = default;
        AlgorithmData(App& app) : AlgorithmData() {
          this->init(app);
        }

        ~AlgorithmData() {
          if(nullptr != initTape) {
            delete initTape;
          }
        }

        IdVector idInY;
        IdVector idInP;
        IdVector idInX;
        IdVector idOutY;
        IdVector idOutP;
        IdVector idOutZ;

        std::vector<RealVector> realCurY;
        std::vector<RealVector> realNextY;

        std::vector<RealVector> realP;
        std::vector<RealVector> realX;

        IdVector idInitX;
        IdVector idInitP;
        Tape* initTape = nullptr;

        void init(App& app) {
          idInY.resize(app.getSizeY());
          idInX.resize(app.getSizeX());
          idOutY.resize(app.getSizeY());
          idOutZ.resize(app.getSizeZ());

          realCurY.resize(app.getNumberOfFunctionals(), RealVector(app.getSizeY()));
          realNextY.resize(app.getNumberOfFunctionals(), RealVector(app.getSizeY()));

          realX.resize(app.getNumberOfFunctionals(), RealVector(app.getSizeX()));

          if(app.getHints() & ApplicationFlags::PIterationIsAvailable) {
            idInP.resize(app.getSizeP());
            idOutP.resize(app.getSizeP());
            realP.resize(app.getNumberOfFunctionals(), RealVector(app.getSizeP()));
          }
        }

        void initInitializationRecording(App& app) {
          idInitX.resize(app.getSizeX());
          idInitP.resize(app.getSizeP());
          initTape = new Tape();
        }
    };

    template<typename T_App>
    struct AlgorithmInterface {
      public:

        using App = CODI_DD(T_App, CODI_T(ApplicationInterface<CODI_ANY>));
        using Type = typename App::Type;

        using Real = RealTraits::Real<Type>;

        using RealVector = std::vector<Real>;

        void run(App& app);

      protected:

        void iterateUntil(App& app, int iteration) {
          while (app.getIteration() < iteration) {
            app.evaluateG();
          }
        }

        struct GetPrimal {
          public:
            RealVector& vec;
            GetPrimal(RealVector& vec) : vec(vec) {}

            void operator()(Type& value, size_t pos) {
              vec[pos] = RealTraits::getValue(value);
            }
        };
    };
  }
}
