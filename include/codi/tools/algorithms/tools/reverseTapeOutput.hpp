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

#include "../../../config.h"
#include "../../../misc/macros.hpp"
#include "../interfaces/algorithmInterface.hpp"
#include "../interfaces/applicationInterface.hpp"
#include "algorithmData.hpp"

/** \copydoc codi::Namespace */
namespace codi {
  namespace algorithms {

    template<typename T_App>
    struct ReverseTapeOutput {
      public:

        using App = CODI_DD(T_App, CODI_T(ApplicationInterface<CODI_ANY>));

        using Type = typename App::Type;
        using Tape = typename Type::Tape;

        using Data = AlgorithmData<App>;
        using RealVector = typename Data::RealVector;
        using IdVector = typename Data::IdVector;

        using ExtFunc = ExternalFunction<Tape>;

      private:

        App& app;
        IdVector ids;
        int iteration;
        OutputType type;
        OutputHints hints;

        ReverseTapeOutput(App& app, IdVector const& ids, int iteration, OutputType type, OutputHints hints) :
          app(app),
          ids(ids),
          iteration(iteration),
          type(type),
          hints(hints) {}

      public:

        static void addReverseOutput(App& app, IdVector const& ids, OutputType type, OutputHints hints) {
          ReverseTapeOutput* out = new ReverseTapeOutput(app, ids, app.getIteration(), type, hints);

          Type::getTape().pushExternalFunction(
                ExtFunc::create(ReverseTapeOutput::reverseOutput, out, ReverseTapeOutput::deleteOutput));
        }

      private:

        static void reverseOutput(Tape* tape, void* data, typename ExtFunc::VectorAccess* adjointInterface) {
          CODI_UNUSED(tape);

          ReverseTapeOutput* out = (ReverseTapeOutput*)data;
          ApplicationIOInterface<Type>* io = out->app.getIOInterface();

          size_t dimMax = adjointInterface->getVectorSize();

          RealVector adj(out->ids.size());
          for(size_t dim = 0; dim < dimMax; dim += 1) {
            for(size_t i = 0; i < out->ids.size(); i += 1) {
              adj[i] = adjointInterface->getAdjoint(out->ids[i], dim);
            }

            io->write(out->type, out->iteration, adj, out->hints, dim);
          }
        }

        static void deleteOutput(Tape* tape, void* data) {
          CODI_UNUSED(tape);

          ReverseTapeOutput* out = (ReverseTapeOutput*)data;
          delete out;
        }
    };
  }
}
