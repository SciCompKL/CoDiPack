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
        int maxIterations; ///< Maximum number of forward iterations.
        std::vector<double> seeding;

        ForwardModeSettings()
            : maxIterations(1000), seeding(1, 1.0)
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
        size_t initSeedingPos;

        ForwardMode() : settings(), initSeedingPos() {}
        ForwardMode(ForwardModeSettings settings) : settings(settings), initSeedingPos() {}

        ForwardModeSettings const* getSettings() const {
          return &settings;
        }

        void run(App& app) {
          ApplicationIOInterface<Type>* io = app.getIOInterface();

          bool isStop = false;
          bool isFinished = false;

          if(ApplicationFlags::InitializationComputesP & app.getHints()) {
            initSeedingPos = 0;
            app.setInitializationHandlingFunction(setGradientInit, this);
          }

          app.initialize();

          if(!(
             1 == settings.seeding.size() ||
             app.getSizeX() * GT::dim == settings.seeding.size())
           ) {
            CODI_EXCEPTION("Seeding of forward mode has the size '%d'. It needs either be one or '%d'.",
                           (int)settings.seeding.size(),
                           (int)app.getSizeX());
          }

          if(!(ApplicationFlags::InitializationComputesP & app.getHints())) {
            app.iterateX(SetGradient(settings.seeding));
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

        static void setGradientInit(Type& value, void* d) {
          ForwardMode* data = (ForwardMode*)d;

          if(data->settings.seeding.size() != 1 && data->initSeedingPos >= data->settings.seeding.size()) {
            CODI_EXCEPTION("Not enough seeding entries.");
          }
          setGradient(value, data->initSeedingPos * GT::dim, data->settings.seeding);
          data->initSeedingPos += 1;
        }

        struct SetGradient {
          public:
            std::vector<double>& seeding;
            SetGradient(std::vector<double>& seeding) : seeding(seeding) {}

            void operator()(Type& value, size_t pos) {
              setGradient(value, pos * GT::dim, seeding);
            }
        };

        static void setGradient(Type& value, size_t pos, std::vector<double>& seeding) {
          for(size_t d = 0; d < GT::dim; d += 1) {
            if(1 == seeding.size()) {
              GT::at(value.gradient(), d) = seeding[0];

            } else {
              GT::at(value.gradient(), d) = seeding[pos + d];
            }
          }
        }
    };
  }
}
