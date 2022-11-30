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

#include "../../../misc/macros.hpp"
#include "../../../config.h"
#include "../interfaces/applicationInterface.hpp"

/** \copydoc codi::Namespace */
namespace codi {
  namespace algorithms {

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

          if(app.getHints() & ApplicationFlags::PStateIsAvailable) {
            idInP.resize(app.getSizeP());
            idOutP.resize(app.getSizeP());
            realP.resize(app.getNumberOfFunctionals(), RealVector(app.getSizeP()));
          }
        }

        void resizeYIn(App& app) {
          if(app.getSizeY() != idInY.size()) {
            idInY.resize(app.getSizeY());

            for(RealVector& vec : realNextY) {
              vec.resize(app.getSizeY());
            }
          }
        }

        void resizeYOut(App& app) {
          if(app.getSizeY() != idOutY.size()) {
            idOutY.resize(app.getSizeY());

            for(RealVector& vec : realCurY) {
              vec.resize(app.getSizeY());
            }
          }
        }

        void initInitializationRecording(App& app) {
          idInitX.resize(app.getSizeX());
          idInitP.resize(app.getSizeP());
          initTape = new Tape();
        }
    };

  }
}
