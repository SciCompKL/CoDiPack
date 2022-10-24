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

    struct ForwardModeSettings {
        int maxIterations;  ///< Maximum number of adjoint iterations.

        ForwardModeSettings()
            : maxIterations(1000)
         {}
    };

    template<typename T_App>
    struct ForwardMode : public AlgorithmInterface<T_App> {
      public:

        using App = CODI_DD(T_App, CODI_T(ApplicationInterface<CODI_ANY>));

        using Type = typename App::Type;

        using Base = AlgorithmInterface<App>;

        using Gradient = typename Type::Gradient;
        using GT = GradientTraits::TraitsImplementation<Gradient>;
        using RealVector = typename Base::RealVector;

        ForwardModeSettings settings;

        ForwardMode() : settings() {}
        ForwardMode(ForwardModeSettings settings) : settings(settings) {}

        ForwardModeSettings const* getSettings() const {
          return &settings;
        }

        void run(App& app) {
          ApplicationIOInterface<Type>* io = app.getIOInterface();

          bool isStop = false;
          bool isFinished = false;

          if(ApplicationFlags::InitializationComputesP & app.getHints()) {
            app.setInitializationHandlingFunction(setGradientInit);
          }

          app.initialize();

          if(!(ApplicationFlags::InitializationComputesP & app.getHints())) {
            app.iterateX(setGradient);
          }

          app.evaluateP();

          app.setInitializationHandlingFunction(nullptr);

          while (!(isFinished || isStop)) {

            app.evaluateG();

            isFinished = app.getIteration() >= settings.maxIterations;
            isStop = app.isStop();
          }

          app.evaluateF();

          std::vector<RealVector> z(GT::dim, RealVector(app.getSizeZ()));
          app.iterateZ(GetGradient(z));

          io->writeZ(app.getIteration(), z, OutputFlags::Derivative | OutputFlags::F | OutputFlags::Final);
        }

        struct GetGradient {
          public:
            std::vector<RealVector>& vec;
            GetGradient(std::vector<RealVector>& vec) : vec(vec) {}

            void operator()(Type& value, size_t pos) {
              for(size_t i = 0; i < GT::dim; i += 1) {
                vec[i][pos] = GT::at(value.getGradient(), i);
              }
            }
        };

        static void setGradientInit(Type& value) {
          value.setGradient(1.0);
        }

        static void setGradient(Type& value, size_t pos) {
          CODI_UNUSED(pos);

          value.setGradient(1.0);
        }
    };
  }
}
