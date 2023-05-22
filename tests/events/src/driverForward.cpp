/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2023 Chair for Scientific Computing (SciComp), University of Kaiserslautern-Landau
 * Homepage: http://www.scicomp.uni-kl.de
 * Contact:  Prof. Nicolas R. Gauger (codi@scicomp.uni-kl.de)
 *
 * Lead developers: Max Sagebaum, Johannes Blühdorn (SciComp, University of Kaiserslautern-Landau)
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
 *  - SciComp, University of Kaiserslautern-Landau:
 *    - Max Sagebaum
 *    - Johannes Blühdorn
 *    - Former members:
 *      - Tim Albring
 */

#include "../include/forwardCallbacks.hpp"
#include "../include/tests/test.hpp"

#ifndef NUMBER
  #error Please define NUMBER as a CoDiPack type.
#endif

int main() {
  using Tape = NUMBER::Tape;
  size_t constexpr dim = codi::GradientTraits::dim<Tape::Gradient>();

  size_t constexpr nInputs = 4;
  size_t constexpr nOutputs = 4;

  auto forwardCallbacks = ForwardCallbacks::registerAll<Tape>();

#ifdef USE_INNER_CALLBACKS
  using InnerTape = Tape::Real::Tape;
  auto innerCallbacks = ForwardCallbacks::registerAll<InnerTape>();
#endif

  NUMBER inputs[nInputs] = {};
  NUMBER outputs[nOutputs] = {};

  size_t constexpr maxRuns = 2;

  for (size_t run = 0; run < maxRuns; run += 1) {
    if (run == maxRuns - 1) { /* last run, deregister all listeners */
      deregisterCallbacks<Tape>(forwardCallbacks);
#ifdef USE_INNER_CALLBACKS
      deregisterCallbacks<InnerTape>(innerCallbacks);
#endif
    }

    std::cout << "# Seed inputs" << std::endl;
    for (size_t i = 0; i < nInputs; ++i) {
      inputs[i] = sin(i + 1);

      for (size_t currentDim = 0; currentDim < dim; ++currentDim) {
        codi::GradientTraits::at(inputs[i].gradient(), currentDim) = cos(i + currentDim * nInputs);
      }

#ifdef USE_INNER_CALLBACKS
      inputs[i].value().setGradient(i + 1);
#endif
    }

    std::cout << "# Run test" << std::endl;
    test<NUMBER>(nInputs, inputs, nOutputs, outputs);
  }

  return 0;
}
