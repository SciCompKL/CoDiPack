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
#include "../../traits/gradientTraits.hpp"
#include "interfaces/algorithmInterface.hpp"

/** \copydoc codi::Namespace */
namespace codi {
  namespace algorithms {

    struct ForwardModeSettings {
        int maxIterations; ///< Maximum number of forward iterations.
        std::vector<double> seeding;

        bool fullJacobian;
        double primalValidationThreshold;

        ForwardModeSettings()
            : maxIterations(1000), seeding(1, 1.0), fullJacobian(false), primalValidationThreshold(1e-10)
         {}
    };

    template<typename T_App>
    struct ForwardMode : public AlgorithmInterface<T_App> {
      public:

        using App = CODI_DD(T_App, CODI_T(ApplicationInterface<CODI_ANY>));

        using Type = typename App::Type;

        using Base = AlgorithmInterface<App>;

        using Real = typename Type::Real;
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
          if(ApplicationHintsFlags::InitializationComputesP & app.getHints() && settings.fullJacobian) {
            CODI_EXCEPTION("Computation of full Jacobian not supported if P can not be recomputed.");
          }

          if(ApplicationHintsFlags::InitializationComputesP & app.getHints()) {
            initSeedingPos = 0;
            app.setInitializationHandlingFunction(setGradientInit, this);
          }

          app.initialize();

          app.setInitializationHandlingFunction(nullptr);

          if(!settings.fullJacobian) {
            runOneTimeMode(app);
          } else {
            runJacobianMode(app);
          }
        }

        void runOneTimeMode(App& app) {
          ApplicationIOInterface<Type>* io = app.getIOInterface();

          if(!(
             1 == settings.seeding.size() ||
             app.getSizeX() * GT::dim == settings.seeding.size())
           ) {
            CODI_EXCEPTION("Seeding of forward mode has the size '%d'. It needs either be one or '%d'.",
                           (int)settings.seeding.size(),
                           (int)app.getSizeX());
          }

          if(!(ApplicationHintsFlags::InitializationComputesP & app.getHints())) {
            app.iterateX(SetGradient(settings.seeding));
          }

          runApp(app);

          std::vector<RealVector> z(GT::dim, RealVector(app.getSizeZ()));
          app.iterateZ(GetGradient(z));

          io->writeZ(app.getIteration(), z, FileOutputHintsFlags::Derivative | FileOutputHintsFlags::F | FileOutputHintsFlags::Final);
        }

        void runJacobianMode(App& app) {
          ApplicationIOInterface<Type>* io = app.getIOInterface();
          CheckpointManagerInterface* cpm = app.getCheckpointInterface();

          FileOutputHints outputHints = FileOutputHintsFlags::Derivative | FileOutputHintsFlags::F | FileOutputHintsFlags::Final | FileOutputHintsFlags::Vector;

          // Create initial checkpoint
          CheckpointHandle* cp = nullptr;
          for(CheckpointHandle* cur : cpm->list()) {
            if(cur->getIteration() == 0) {
              cp = cur;
              break;
            }
          }
          if(nullptr == cp) {
            cp = cpm->create();
          }


          RealVector zValue(app.getSizeZ());
          std::vector<RealVector> zGrad(GT::dim, RealVector(app.getSizeZ()));

          size_t sizeX = app.getSizeX();
          for(size_t curX = 0; curX < sizeX; curX += GT::dim) {
            app.print(StringUtil::format("Computing %d/%d (Vec: %d)\n", (int)(curX + 1), (int)sizeX, (int)GT::dim));
            app.iterateX(SetGradientAtPos(1.0, curX, GT::dim));

            runApp(app);

            if(curX + GT::dim >= sizeX) {
                  zGrad.resize(sizeX - curX); // Resize to last vector dimension
            }
            app.iterateZ(GetGradient(zGrad));
            if(0 == curX) {
              app.iterateZ(GetValue(zValue));
            } else {
              size_t errors = 0;
              app.iterateZ(ValidateValue(zValue, errors, settings.primalValidationThreshold));
              if(0 != errors) {
                app.print(StringUtil::format("Warning: Primal changed in '%d' places in the '%d' run.\n", (int)errors, (int)curX));
              }
            }

            io->writeZ(app.getIteration(), zGrad, outputHints, curX);

            cpm->load(cp);
          }
        }

        void runApp(App& app) {
          app.evaluateP();

          bool isStop = false;
          bool isFinished = false;
          while (!(isFinished || isStop)) {

            app.evaluateG();

            isFinished = app.getIteration() >= settings.maxIterations;
            isStop = app.isStop();
          }

          app.evaluateF();
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

        struct GetValue {
          public:
            RealVector& vec;
            GetValue(RealVector& vec) : vec(vec) {}

            void operator()(Type& value, size_t pos) {
              vec[pos] = value.getValue();
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

        struct SetGradientAtPos {
          public:
            Real gradient;
            size_t pos;
            size_t vecDim;
            SetGradientAtPos(Real const& gradient, size_t const& pos, size_t const vecDim) : gradient(gradient), pos(pos), vecDim(vecDim) {}

            void operator()(Type& value, size_t pos) {
              if(this->pos <= pos && pos < this->pos + vecDim) {
                size_t dim = pos - this->pos;
                GT::at(value.gradient(), dim) = gradient;
              } else {
                value.setGradient(Gradient());
              }
            }
        };

        struct ValidateValue {
          public:
            RealVector const& vec;
            size_t& errors;
            double threshold;
            ValidateValue(RealVector const& vec, size_t errors, double threshold) :
              vec(vec), errors(errors), threshold(threshold) {}

            void operator()(Type& value, size_t pos) {
              Real diff = value.getValue() - vec[pos];
              if(Real() != vec[pos]) {
                diff = diff / vec[pos];
              }
              if(abs(diff) >= threshold) {
                errors += 1;
              }
            }
        };

    };
  }
}
