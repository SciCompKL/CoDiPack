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

#ifdef SECOND_ORDER
      using InnerTape = typename Tape::Real::Tape;
      auto innerCallbacks = ForwardCallbacks::registerAll<InnerTape>();
#endif

      Number inputs[nInputs] = {};
      Number outputs[nOutputs] = {};

      size_t constexpr maxRuns = 3;

      for (size_t run = 0; run < maxRuns; run += 1) {
        if (run == maxRuns - 1) { /* last run, deregister all listeners */
          deregisterCallbacks<Tape>(reverseCallbacks);
#ifdef SECOND_ORDER
          deregisterCallbacks<InnerTape>(innerCallbacks);
#endif
        }

        tape.reset();

        tape.setActive();

        std::cout << "# Register inputs" << std::endl;
        for (size_t i = 0; i < nInputs; ++i) {
          inputs[i] = sin(i + 1);

#ifdef SECOND_ORDER
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
#ifdef SECOND_ORDER
      ForwardCallbacks::registerAll<InnerTape>();
#endif

      tape.resetHard();
    }

    virtual void evaluate(Tape& tape) {
      tape.evaluate();
    }
};
