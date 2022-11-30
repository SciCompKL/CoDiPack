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

#include <algorithm>
#include <vector>

#include "../../config.h"
#include "../../expressions/lhsExpressionInterface.hpp"
#include "../../misc/macros.hpp"
#include "../../misc/stringUtil.hpp"
#include "interfaces/algorithmInterface.hpp"

/** \copydoc codi::Namespace */
namespace codi {
  namespace algorithms {

    struct CheckpointTestSettings {
      public:
        int checkpointIter;
        std::vector<int> compareIter;

        double maxRelativeError;

        bool forceWrite;

        CheckpointTestSettings() : checkpointIter(10), compareIter(), maxRelativeError(1e-12), forceWrite(false) {
          compareIter.push_back(10);
          compareIter.push_back(20);
          compareIter.push_back(30);
        }
    };

    template<typename T_App>
    struct CheckpointTest : public AlgorithmInterface<T_App> {
      public:

        using App = CODI_DD(T_App, CODI_T(ApplicationInterface<CODI_ANY>));

        using Type = typename App::Type;
        using Real = typename App::Real;

        using Base = AlgorithmInterface<App>;
        using Res = typename App::Res;
        using RealVector = typename Base::RealVector;

        CheckpointTestSettings settings;

        CheckpointTest(CheckpointTestSettings settings) : settings(settings) {}

        CheckpointTestSettings const* getSettings() const {
          return &settings;
        }

        void run(App& app) {
          CheckpointManagerInterface* cm = app.getCheckpointInterface();
          ApplicationIOInterface<Type>* io = app.getIOInterface();

          if(app.getHints() & ApplicationHintsFlags::InitializationRequired) {
            app.initialize();
          }

          // Pare settings
          std::sort(settings.compareIter.begin(), settings.compareIter.end());
          settings.compareIter.erase(std::unique(settings.compareIter.begin(), settings.compareIter.end()),
                                     settings.compareIter.end());
          validateSettings(app);

          if (app.getIteration() < settings.checkpointIter) {
            app.print(StringUtil::format("Iterating to checkpoint iteration %d.\n", settings.checkpointIter));
            Base::iterateUntil(app, settings.checkpointIter);
          }

          app.print("Creating checkpoint.\n");
          CheckpointHandle* cp = cm->create();

          int nComparisons = settings.compareIter.size();
          std::vector<RealVector> vectors(nComparisons);

          for (size_t i = 0; i < settings.compareIter.size(); i += 1) {
            int compareIter = settings.compareIter[i];

            if (app.getIteration() < compareIter) {
              app.print(StringUtil::format("Iterating to comparison iteration %d.\n", compareIter));
              Base::iterateUntil(app, compareIter);
            }

            app.print(StringUtil::format("Getting solution at iteration %d.\n", compareIter));
            vectors[i].resize(app.getSizeY());
            app.iterateY(typename Base::GetPrimal(vectors[i]));
          }

          app.print(StringUtil::format("Restoring checkpoint at %d.\n", settings.checkpointIter));
          cm->load(cp);
          bool correctIteration = settings.checkpointIter == app.getIteration();
          app.print(StringUtil::format("Iteration is correctly reset %d. (%d == %d).\n", (int)correctIteration,
                                       settings.checkpointIter, app.getIteration()));

          RealVector curPrimal;
          for (size_t i = 0; i < settings.compareIter.size(); i += 1) {
            int compareIter = settings.compareIter[i];

            if (app.getIteration() < compareIter) {
              app.print(StringUtil::format("Iterating to comparison iteration %d.\n", compareIter));
              Base::iterateUntil(app, compareIter);
            }

            app.print(StringUtil::format("Getting solution at iteration %d.\n", compareIter));
            curPrimal.resize(app.getSizeY());
            app.iterateY(typename Base::GetPrimal(curPrimal));

            app.print(StringUtil::format("Comparing current solution with stored one ..", compareIter));
            double largestError;
            int errorCount;
            compareVectors(vectors[i], curPrimal, largestError, errorCount);
            if(0 == errorCount && !settings.forceWrite) {
              app.print("OK\n");
            } else {
              app.print(StringUtil::format("found %d errors, largest is %0.6e.\n", errorCount, largestError));
              app.print(StringUtil::format("Writing vectors."));
              io->writeY(compareIter, vectors[i],
                         FileOutputHintsFlags::Primal | FileOutputHintsFlags::Intermediate | FileOutputHintsFlags::G | FileOutputHintsFlags::V1);
              io->writeY(compareIter, curPrimal,
                         FileOutputHintsFlags::Primal | FileOutputHintsFlags::Intermediate | FileOutputHintsFlags::G | FileOutputHintsFlags::V2);
            }
          }
        }

      protected:

        inline void compareVectors(RealVector const& v1, RealVector const& v2, double& largestError, int& errorCount) {
          largestError = -1e300;
          errorCount = 0;

          if (v1.size() == v2.size()) {
            for (size_t pos = 0; pos < v2.size(); pos += 1) {
              double base = RealTraits::getPassiveValue(v1[pos]);
              double diff = abs(RealTraits::getPassiveValue(v2[pos]) - base);

              double rel = diff / abs(base);
              if (rel > settings.maxRelativeError) {
                errorCount += 1;
                largestError = max(largestError, rel);
              }
            }
          } else {
            errorCount = max(v1.size(), v2.size());
            largestError = 1.0;
          }
        }

        inline void validateSettings(App& app) {
          int curIteration = app.getIteration();

          bool hasError = false;
          if (curIteration > settings.checkpointIter) {
            app.print(StringUtil::format("Error: Current app iteration(%d) is behind the checkpoint iteration(%d).",
                                         curIteration, settings.checkpointIter));
            hasError = true;
          }

          for (int compareIter : settings.compareIter) {
            if (compareIter < settings.checkpointIter) {
              app.print(
                  StringUtil::format("Error: Iteration for comparison(%d) is before the checkpoint iteration(%d).",
                                     compareIter, settings.checkpointIter));
              hasError = true;
            }
          }

          if (hasError) {
            exit(-1);
          }
        }
    };
  }
}
