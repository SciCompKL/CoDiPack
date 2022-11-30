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
#include "../../traits/realTraits.hpp"
#include "interfaces/algorithmInterface.hpp"

/** \copydoc codi::Namespace */
namespace codi {
  namespace algorithms {

    struct FiniteDifferenceEvaluationSettings {
        int maxIterations;  ///< Maximum number of adjoint iterations.

        bool fullJacobian;

        std::vector<double> stepSizes;
        bool relativeStepSize;
        bool validateBase;
        double primalValidationThreshold;

        bool writePrimal;

        FiniteDifferenceEvaluationSettings()
            : maxIterations(1000), fullJacobian(false), stepSizes(1),
              relativeStepSize(true), validateBase(true), primalValidationThreshold(1e-10), writePrimal(false)
         {
          stepSizes[0] = 0.1; // 10% distortion
        }
    };

    template<typename T_App>
    struct FiniteDifferenceEvaluation : public AlgorithmInterface<T_App> {
      public:

        using App = CODI_DD(T_App, CODI_T(ApplicationInterface<CODI_ANY>));

        using Type = typename App::Type;

        using Base = AlgorithmInterface<App>;

        using Real = RealTraits::Real<Type>;
        using RealVector = typename Base::RealVector;

        FiniteDifferenceEvaluationSettings settings;

        FiniteDifferenceEvaluation() : settings() {}
        FiniteDifferenceEvaluation(FiniteDifferenceEvaluationSettings settings) : settings(settings) {}

        FiniteDifferenceEvaluationSettings const* getSettings() const {
          return &settings;
        }

        void run(App& app) {
          if(ApplicationHintsFlags::InitializationComputesP & app.getHints() && settings.fullJacobian) {
            CODI_EXCEPTION("Computation of full Jacobian not supported if P can not be recomputed.");
          }

          if(ApplicationHintsFlags::InitializationComputesP & app.getHints()) {
            // TODO: Implement direct init.
          }

          app.initialize();

          if(!settings.fullJacobian) {
            runOneTimeMode(app);
          } else {
            runJacobianMode(app);
          }
        }

        void runOneTimeMode(App& app) {
          CODI_EXCEPTION("Not implemented."); // TODO: implement
        }

        void runJacobianMode(App& app) {
          ApplicationIOInterface<Type>* io = app.getIOInterface();
          CheckpointManagerInterface* cpm = app.getCheckpointInterface();

          FileOutputHints primalHints = FileOutputHintsFlags::Primal | FileOutputHintsFlags::F | FileOutputHintsFlags::Final;
          FileOutputHints gradientHints = FileOutputHintsFlags::Derivative | FileOutputHintsFlags::F | FileOutputHintsFlags::Final | FileOutputHintsFlags::Vector;

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

          RealVector xBase(app.getSizeX());
          RealVector zBase(app.getSizeZ());
          RealVector zGrad(app.getSizeZ());

          app.iterateX(GetValue(xBase));

          app.print(StringUtil::format("Computing base.\n"));
          runApp(app);
          app.iterateZ(GetValue(zBase));
          if(settings.writePrimal) {
            io->changeFolder("base");
            io->writeZ(app.getIteration(), zBase, primalHints);
          }
          cpm->load(cp);

          // Validate base
          if(settings.validateBase) {
            app.print(StringUtil::format("Validating base.\n"));
            runApp(app);
            size_t errors = 0;
            app.iterateZ(ValidateValue(zBase, errors, settings.primalValidationThreshold));
            if(0 != errors) {
              io->writeZ(app.getIteration(), zBase,  primalHints | FileOutputHintsFlags::V1);
              app.iterateZ(GetValue(zBase));
              io->writeZ(app.getIteration(), zBase,  primalHints | FileOutputHintsFlags::V2);
              CODI_EXCEPTION("Error: Primal changed in '%d' places.\n", (int)errors);
            }
            cpm->load(cp);
          }

          size_t sizeX = app.getSizeX();
          for(size_t curStep = 0; curStep < settings.stepSizes.size(); curStep += 1) {
            io->changeFolder(StringUtil::format("step_%04d", (int)curStep));

            for(size_t curX = 0; curX < sizeX; curX += 1) {
              app.print(StringUtil::format("Computing step: %d/%d (%.6e) input: %d/%d.\n",
                                           (int)(curStep + 1), (int) settings.stepSizes.size(),
                                           settings.stepSizes[curStep],
                                           (int)(curX + 1), (int)sizeX));
              Real actualStepSize;
              app.iterateX(SetPrimalOffsetAtPos(xBase, curX, settings.stepSizes[curStep], settings.relativeStepSize, actualStepSize));

              runApp(app);
              app.iterateZ(GetValue(zGrad));

              if(settings.writePrimal) {
                io->writeZ(app.getIteration(), zGrad, primalHints, curX);
              }

              computeGrad(zGrad, zBase, actualStepSize, settings.relativeStepSize);
              io->writeZ(app.getIteration(), zGrad, gradientHints, curX);

              cpm->load(cp);
            }
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

        void computeGrad(RealVector& grad, RealVector const& base, Real const& step, bool relative) {
          for(size_t pos = 0; pos < grad.size(); pos += 1) {
            Real diff = grad[pos] - base[pos];
            grad[pos] = diff / step;
          }
        }

        struct GetValue {
          public:
            RealVector& vec;
            GetValue(RealVector& vec) : vec(vec) {}

            void operator()(Type& value, size_t pos) {
              vec[pos] = RealTraits::getValue(value);
            }
        };

        struct SetPrimalOffsetAtPos {
          public:
            RealVector const& base;
            size_t pos;
            Real stepSize;
            bool relative;
            Real& actualStepSize;
            SetPrimalOffsetAtPos(RealVector const& base, size_t const& pos, Real const& stepSize, bool const& relative,
                                 Real& actualStepSize) :
              base(base), pos(pos), stepSize(stepSize), relative(relative), actualStepSize(actualStepSize) {}

            void operator()(Type& value, size_t pos) {
              value = base[pos];
              if(this->pos == pos) {
                actualStepSize = stepSize;
                if(relative && Real() != RealTraits::getValue(value)) {
                  actualStepSize *= RealTraits::getValue(value);
                }
                value += actualStepSize;
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
              Real diff = RealTraits::getValue(value) - vec[pos];
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
