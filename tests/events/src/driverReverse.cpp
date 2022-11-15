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

#include "../include/forwardCallbacks.hpp"
#include "../include/reverseCallbacks.hpp"
#include "../include/tests/test.hpp"

#ifndef NUMBER
  #error Please define NUMBER as a CoDiPack type.
#endif

int main() {
  using Tape = NUMBER::Tape;
  size_t constexpr dim = codi::GradientTraits::dim<Tape::Gradient>();

  auto& tape = NUMBER::getTape();

  size_t constexpr nInputs = 4;
  size_t constexpr nOutputs = 4;

  ReverseCallbacks::registerAll<Tape>();

#ifdef SECOND_ORDER
  using InnerTape = Tape::Real::Tape;
  ForwardCallbacks::registerAll<InnerTape>();
#endif

  NUMBER inputs[nInputs] = {};
  NUMBER outputs[nOutputs] = {};

  size_t constexpr maxRuns = 2;

  for (size_t run = 0; run < maxRuns; run += 1) {
    tape.reset();

    tape.setActive();

    for (size_t i = 0; i < nInputs; ++i) {
      inputs[i] = sin(i + 1);

#ifdef SECOND_ORDER
      inputs[i].value().setGradient(i + 1);
#endif

      tape.registerInput(inputs[i]);
    }

    test<NUMBER>(nInputs, inputs, nOutputs, outputs);

    for (size_t j = 0; j < nOutputs; ++j) {
      tape.registerOutput(outputs[j]);
    }

    tape.setPassive();

    for (size_t j = 0; j < nOutputs; ++j) {
      for (size_t currentDim = 0; currentDim < dim; ++currentDim) {
        codi::GradientTraits::at(outputs[j].gradient(), currentDim) = cos(dim * j + currentDim);
      }
    }

    tape.evaluate();

    ReverseCallbacks::GlobalStatementCounters<Tape>::assertEqual();
  }

  tape.resetHard();

  return 0;
}
