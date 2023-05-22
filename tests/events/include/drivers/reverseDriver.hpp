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
#pragma once

#include "../forwardCallbacks.hpp"
#include "../reverseCallbacks.hpp"
#include "../tests/test.hpp"

template<typename Number>
struct ReverseDriver {
    using Tape = typename Number::Tape;

    void run() {
      size_t constexpr dim = codi::GradientTraits::dim<typename Tape::Gradient>();

      auto& tape = Number::getTape();

      size_t constexpr nInputs = 4;
      size_t constexpr nOutputs = 4;

      auto reverseCallbacks = ReverseCallbacks::registerAll<Tape>();

#ifdef USE_INNER_CALLBACKS
      using InnerTape = typename Tape::Real::Tape;
      auto innerCallbacks = ForwardCallbacks::registerAll<InnerTape>();
#endif

      Number inputs[nInputs] = {};
      Number outputs[nOutputs] = {};

      size_t constexpr maxRuns = 3;

      for (size_t run = 0; run < maxRuns; run += 1) {
        if (run == maxRuns - 1) { /* last run, deregister all listeners */
          deregisterCallbacks<Tape>(reverseCallbacks);
#ifdef USE_INNER_CALLBACKS
          deregisterCallbacks<InnerTape>(innerCallbacks);
#endif
        }

        tape.reset();

        tape.setActive();

        std::cout << "# Register inputs" << std::endl;
        for (size_t i = 0; i < nInputs; ++i) {
          inputs[i] = sin(i + 1);

#ifdef USE_INNER_CALLBACKS
          inputs[i].value().setGradient(i + 1);
#endif

          tape.registerInput(inputs[i]);
        }

        std::cout << "# Run test" << std::endl;
        test<Number>(nInputs, inputs, nOutputs, outputs);

        std::cout << "# Register outputs" << std::endl;
        for (size_t j = 0; j < nOutputs; ++j) {
          tape.registerOutput(outputs[j]);
        }

        tape.setPassive();

        for (size_t j = 0; j < nOutputs; ++j) {
          for (size_t currentDim = 0; currentDim < dim; ++currentDim) {
            codi::GradientTraits::at(outputs[j].gradient(), currentDim) = cos(j + currentDim * nOutputs);
          }
        }

        std::cout << "# Tape evaluate" << std::endl;
        evaluate(tape);

        ReverseCallbacks::GlobalStatementCounters<Tape>::assertEqual();
      }

      /* re-register for testing resetHard */
      ReverseCallbacks::registerAll<Tape>();
#ifdef USE_INNER_CALLBACKS
      ForwardCallbacks::registerAll<InnerTape>();
#endif

      tape.resetHard();
    }

    virtual void evaluate(Tape& tape) {
      tape.evaluate();
    }
};
