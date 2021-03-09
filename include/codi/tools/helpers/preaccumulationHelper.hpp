
#pragma once

#include <vector>

#include "../../aux/exceptions.hpp"
#include "../../config.h"
#include "../../expressions/lhsExpressionInterface.hpp"
#include "../../tapes/interfaces/fullTapeInterface.hpp"
#include "../../traits/gradientTraits.hpp"
#include "../../traits/tapeTraits.hpp"
#include "../algorithms.hpp"
#include "../data/jacobian.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  template<typename _Type, typename = void>
  struct PreaccumulationHelper {
    public:

      using Type = CODI_DECLARE_DEFAULT(_Type,
                                        CODI_TEMPLATE(LhsExpressionInterface<double, double, CODI_ANY, CODI_ANY>));

      using Real = typename Type::Real;
      using Identifier = typename Type::Identifier;
      using Gradient = typename Type::Gradient;

      using Tape = CODI_DECLARE_DEFAULT(typename Type::Tape,
                                        CODI_TEMPLATE(FullTapeInterface<double, double, int, CODI_ANY>));
      using Position = typename Tape::Position;

      std::vector<Identifier> inputData;
      std::vector<Identifier> outputData;
      std::vector<Type*> outputValues;
      Position startPos;

      std::vector<Gradient> storedAdjoints;

      JacobianCountNonZerosRow<Real> jacobie;

      PreaccumulationHelper()
          : inputData(), outputData(), outputValues(), startPos(), storedAdjoints(), jacobie(0, 0) {}

      template<typename... Inputs>
      void addInput(Inputs const&... inputs) {
        Tape& tape = Type::getGlobalTape();

        if (tape.isActive()) {
          addInputRec(inputs...);
        }
      }

      template<typename... Inputs>
      void start(Inputs const&... inputs) {
        Tape& tape = Type::getGlobalTape();

        if (tape.isActive()) {
          inputData.clear();
          outputData.clear();
          outputValues.clear();

          startPos = tape.getPosition();

          addInputRec(inputs...);
        }
      }

      template<typename... Outputs>
      void addOutput(Outputs&... outputs) {
        Tape& tape = Type::getGlobalTape();

        if (tape.isActive()) {
          addOutputRec(outputs...);
        }
      }

      template<typename... Outputs>
      void finish(bool const storeAdjoints, Outputs&... outputs) {
        Tape& tape = Type::getGlobalTape();

        if (tape.isActive()) {
          addOutputRec(outputs...);

          if (storeAdjoints) {
            storeInputAdjoints();
          }

          doPreaccumulation();

          if (storeAdjoints) {
            restoreInputAdjoints();
          }
        }
      }

    private:

      void addInputLogic(Type const& input) {
        Identifier const& identifier = input.getIdentifier();
        if (0 != identifier) {
          inputData.push_back(identifier);
        }
      }

      /// Terminator for the recursive implementation
      void addInputRec() {
        // terminator implementation
      }

      template<typename... Inputs>
      void addInputRec(Type const& input, Inputs const&... r) {
        addInputLogic(input);
        addInputRec(r...);
      }

      void addOutputLogic(Type& output) {
        Identifier const& identifier = output.getIdentifier();
        if (0 != identifier) {
          outputData.push_back(identifier);
          outputValues.push_back(&output);
        }
      }

      /// Terminator for the recursive implementation
      void addOutputRec() {
        // terminator implementation
      }

      template<typename... Outputs>
      void addOutputRec(Type& output, Outputs&... r) {
        addOutputLogic(output);
        addOutputRec(r...);
      }

      void storeInputAdjoints() {
        Tape& tape = Type::getGlobalTape();

        if (storedAdjoints.size() < inputData.size()) {
          storedAdjoints.resize(inputData.size());
        }

        for (size_t i = 0; i < inputData.size(); ++i) {
          Identifier index = inputData[i];
          Gradient& adjoint = tape.gradient(index);
          storedAdjoints[i] = adjoint;
          adjoint = Gradient();
        }
      }

      void restoreInputAdjoints() {
        Tape& tape = Type::getGlobalTape();

        for (size_t i = 0; i < inputData.size(); ++i) {
          Identifier index = inputData[i];
          tape.gradient(index) = storedAdjoints[i];
        }
      }

      void doPreaccumulation() {
        // perform the accumulation of the tape part
        Tape& tape = Type::getGlobalTape();

        Position endPos = tape.getPosition();
        if (jacobie.getM() != outputData.size() || jacobie.getN() != inputData.size()) {
          jacobie.resize(outputData.size(), inputData.size());
        }

        Algorithms<Type, false>::computeJacobian(startPos, endPos, inputData.data(), inputData.size(),
                                                 outputData.data(), outputData.size(), jacobie);

        // store the Jacobian matrix
        tape.resetTo(startPos);

        for (size_t curOut = 0; curOut < outputData.size(); ++curOut) {
          Type& value = *outputValues[curOut];
          if (0 != jacobie.nonZerosRow(curOut)) {
            int nonZerosLeft = jacobie.nonZerosRow(curOut);
            jacobie.nonZerosRow(curOut) = 0;

            // we need to use here the value of the gradient data such that it is correctly deleted.
            Identifier lastIdentifier = value.getIdentifier();
            bool staggeringActive = false;
            int curIn = 0;

            // push statements as long as there are non zeros left
            // if there are more than MaxStatementIntValue non zeros, then we need to stagger the
            // statement pushes
            while (nonZerosLeft > 0) {
              // calculate the number of Jacobians for this statement
              int jacobiesForStatement = nonZerosLeft;
              if (jacobiesForStatement > (int)Config::MaxArgumentSize) {
                jacobiesForStatement = (int)Config::MaxArgumentSize - 1;
                if (staggeringActive) {  // Space is used up but we need one Jacobian for the staggering
                  jacobiesForStatement -= 1;
                }
              }
              nonZerosLeft -= jacobiesForStatement;  // Update non zeros so that we know if it is the last round

              Identifier storedIdentifier = lastIdentifier;
              tape.storeManual(value.getValue(), lastIdentifier, jacobiesForStatement + (int)staggeringActive);
              if (staggeringActive) {  // Not the first staggering so push the last output
                tape.pushJacobiManual(1.0, 0.0, storedIdentifier);
              }

              // push the rest of the Jacobians for the statement
              while (jacobiesForStatement > 0) {
                if (Real() != (Real)jacobie(curOut, curIn)) {
                  tape.pushJacobiManual(jacobie(curOut, curIn), 0.0, inputData[curIn]);
                  jacobiesForStatement -= 1;
                }
                curIn += 1;
              }

              staggeringActive = true;
            }

            value.getIdentifier() = lastIdentifier; /* now set gradient data for the real output value */
          } else {
            // disable tape index since there is no dependency
            tape.destroyIdentifier(value.value(), value.getIdentifier());
          }
        }
      }
  };

  /**
   * @brief Helper implementation of the same interface as the PreaccumulationHelper for forward AD tapes.
   *
   * This implementation does nothing in all methods.
   */
  struct PreaccumulationHelperNoOpBase {
    public:

      /// Does nothing
      template<typename... Inputs>
      void addInput(Inputs const&... inputs) {
        CODI_UNUSED(inputs...);
        // do nothing
      }

      /// Does nothing
      template<typename... Inputs>
      void start(Inputs const&... inputs) {
        CODI_UNUSED(inputs...);
        // do nothing
      }

      /// Does nothing
      template<typename... Outputs>
      void addOutput(Outputs&... outputs) {
        CODI_UNUSED(outputs...);
        // do nothing
      }

      /// Does nothing
      template<typename... Outputs>
      void finish(bool const storeAdjoints, Outputs&... outputs) {
        CODI_UNUSED(storeAdjoints, outputs...);
        // do nothing
      }
  };

#ifndef DOXYGEN_DISABLE
  /// Spcialize PreaccumulationHelper for forward tapes
  template<typename Type>
  struct PreaccumulationHelper<Type, TapeTraits::EnableIfForwardTape<typename Type::Tape>>
      : public PreaccumulationHelperNoOpBase {};

  /// Spcialize PreaccumulationHelper for doubles
  template<>
  struct PreaccumulationHelper<double, void> : public PreaccumulationHelperNoOpBase {};
#endif
}
