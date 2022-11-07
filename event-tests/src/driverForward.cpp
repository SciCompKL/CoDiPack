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
#include "../include/test.hpp"

#ifndef NUMBER
  #error Please define NUMBER as a CoDiPack type.
#endif

int main() {

  using Tape = NUMBER::Tape;
  size_t constexpr dim = codi::GradientTraits::dim<Tape::Gradient>();

  size_t constexpr nInputs = 4;
  size_t constexpr nOutputs = 4;

  registerForwardCallbacks<Tape>();

  NUMBER inputs[nInputs] = {};
  NUMBER outputs[nOutputs] = {};

  for (int offset = 0; offset < nInputs; offset += dim) {

    for (size_t i = 0; i < nInputs; ++i) {
      inputs[i] = sin(i + 1);
    }

    for (size_t currentDim = 0; currentDim < dim && offset + currentDim < nInputs; ++currentDim) {
      codi::GradientTraits::at(inputs[offset + currentDim].gradient(), currentDim) = 1.0;
    }

    test<NUMBER>(nInputs, inputs, nOutputs, outputs);
  }

  return 0;
}
