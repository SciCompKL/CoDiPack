
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

  /**
   * @brief Stores the Jacobian matrix for a code section.
   *
   * The preaccumulation of a code section corresponds to the process of replacing the recorded tape entries with the
   * Jacobian matrix of that section. If the code part is defined by the function \f$ f \f$, then the Jacobian
   * \f$ \frac{df}{dx} \f$ is computed by the Preaccumulation helper and stored on the tape.
   *
   * The preaccumulation of a code part is beneficial if it is complicated to compute but has only a few inputs and
   * outputs. If the computation requires 200 statements with a total of 600 arguments, the storage for this would be
   * 7400 byte on a Jacobian tape. If the function has only two input arguments and two output arguments, the storage
   * for the Jacobian matrix of this function would require 50 byte.
   *
   * The procedure for the preaccumulation of a code section is as follows.
   *
   * \snippet documentation/examples/Example_15_Preaccumulation_of_code_parts.cpp Preaccumulation region
   *
   * The preaccumulation helper can be used multiple times, the start routine resets the state such that multiple
   * evaluations are possible. This improves the performance of the helper since stack allocations are only performed
   * once.
   *
   * @tparam T_Type  The CoDiPack type on which the evaluations take place.
   */
  template<typename T_Type, typename = void>
  struct PreaccumulationHelper {
    public:

      /// See PreaccumulationHelper.
      using Type = CODI_DECLARE_DEFAULT(T_Type,
                                        CODI_TEMPLATE(LhsExpressionInterface<double, double, CODI_ANY, CODI_ANY>));

      using Real = typename Type::Real;              ///< See LhsExpressionInterface.
      using Identifier = typename Type::Identifier;  ///< See LhsExpressionInterface.
      using Gradient = typename Type::Gradient;      ///< See LhsExpressionInterface.

      /// See LhsExpressionInterface.
      using Tape = CODI_DECLARE_DEFAULT(typename Type::Tape,
                                        CODI_TEMPLATE(FullTapeInterface<double, double, int, CODI_ANY>));
      using Position = typename Tape::Position;  ///< See PositionalEvaluationTapeInterface.

      std::vector<Identifier> inputData;   ///< List of input identifiers. Can be added manually after start() was
                                           ///< called.
      std::vector<Identifier> outputData;  ///< List of output identifiers. Can be added manually before finish() is
                                           ///< called. Has to be in sync with outputValues.
      std::vector<Type*> outputValues;     ///< List of output value pointers. Can be added manually before finish() is
                                           ///< called. Has to be in sync with outputData.

    protected:

      Position startPos;                       ///< Starting position for the region.
      std::vector<Gradient> storedAdjoints;    ///< If adjoints of inputs should be stored, before the preaccumulation.
      JacobianCountNonZerosRow<Real> jacobiean;  ///< Jacobian for the preaccumulation.

    public:

      /// Constructor
      PreaccumulationHelper()
          : inputData(), outputData(), outputValues(), startPos(), storedAdjoints(), jacobiean(0, 0) {}

      /// Add multiple additional inputs. Inputs need to be of type `Type`. Called after start().
      template<typename... Inputs>
      void addInput(Inputs const&... inputs) {
        Tape& tape = Type::getTape();

        if (tape.isActive()) {
          addInputRecursive(inputs...);
        }
      }

      /// Starts a preaccumulation region. Resets the internal state. See `addInputs()` for inputs.
      template<typename... Inputs>
      void start(Inputs const&... inputs) {
        Tape& tape = Type::getTape();

        if (tape.isActive()) {
          inputData.clear();
          outputData.clear();
          outputValues.clear();

          startPos = tape.getPosition();

          addInputRecursive(inputs...);
        }
      }

      /// Add multiple additional outputs. Outputs need to be of type `Type`. Called before finish().
      template<typename... Outputs>
      void addOutput(Outputs&... outputs) {
        Tape& tape = Type::getTape();

        if (tape.isActive()) {
          addOutputRecursive(outputs...);
        }
      }

      /// Finish the preaccumulation region and perform the preaccumulation. See `addOutput()` for outputs.
      template<typename... Outputs>
      void finish(bool const storeAdjoints, Outputs&... outputs) {
        Tape& tape = Type::getTape();

        if (tape.isActive()) {
          addOutputRecursive(outputs...);

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

      /// Terminator for the recursive implementation.
      void addInputRecursive() {
        // Terminator implementation.
      }

      template<typename... Inputs>
      void addInputRecursive(Type const& input, Inputs const&... r) {
        addInputLogic(input);
        addInputRecursive(r...);
      }

      void addOutputLogic(Type& output) {
        Identifier const& identifier = output.getIdentifier();
        if (0 != identifier) {
          outputData.push_back(identifier);
          outputValues.push_back(&output);
        }
      }

      /// Terminator for the recursive implementation.
      void addOutputRecursive() {
        // Terminator implementation.
      }

      template<typename... Outputs>
      void addOutputRecursive(Type& output, Outputs&... r) {
        addOutputLogic(output);
        addOutputRecursive(r...);
      }

      void storeInputAdjoints() {
        Tape& tape = Type::getTape();

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
        Tape& tape = Type::getTape();

        for (size_t i = 0; i < inputData.size(); ++i) {
          Identifier index = inputData[i];
          tape.gradient(index) = storedAdjoints[i];
        }
      }

      void doPreaccumulation() {
        // Perform the accumulation of the tape part.
        Tape& tape = Type::getTape();

        Position endPos = tape.getPosition();
        if (jacobiean.getM() != outputData.size() || jacobiean.getN() != inputData.size()) {
          jacobiean.resize(outputData.size(), inputData.size());
        }

        Algorithms<Type, false>::computeJacobian(startPos, endPos, inputData.data(), inputData.size(),
                                                 outputData.data(), outputData.size(), jacobiean);

        // Store the Jacobian matrix.
        tape.resetTo(startPos);

        for (size_t curOut = 0; curOut < outputData.size(); ++curOut) {
          Type& value = *outputValues[curOut];
          if (0 != jacobiean.nonZerosRow(curOut)) {
            int nonZerosLeft = jacobiean.nonZerosRow(curOut);
            jacobiean.nonZerosRow(curOut) = 0;

            // We need to initialize with the output's current identifier such that it is correctly deleted in
            // storeManual.
            Identifier lastIdentifier = value.getIdentifier();
            bool staggeringActive = false;
            int curIn = 0;

            // Push statements as long as there are nonzeros left.
            // If there are more than MaxStatementIntValue nonzeros, then we need to stagger the
            // statement pushes:
            // e.g. The reverse mode of w = f(u0, ..., u530) which is \bar u_i += df/du_i * \bar w for i = 0 ... 530 is
            //      separated into
            //        Statement 1:
            //          \bar u_i += df/du_i * \bar t_1 for i = 0 ... 253   (254 entries)
            //        Statement 2:
            //          \bar t_1 += \bar w                                 (1 entry)
            //          \bar u_i += df/du_i * \bar t_2 for i = 254 ... 506 (253 entries)
            //        Statement 3:
            //          \bar t_2 += \bar w                                 (1 entry)
            //          \bar u_i += df/du_i * \bar w for i = 507 ... 530   (24 entries)
            //
            while (nonZerosLeft > 0) {
              // Calculate the number of Jacobians for this statement.
              int jacobiansForStatement = nonZerosLeft;
              if (jacobiansForStatement > (int)Config::MaxArgumentSize) {
                jacobiansForStatement = (int)Config::MaxArgumentSize - 1;
                if (staggeringActive) {  // Except in the first round, one Jacobian is reserved for the staggering.
                  jacobiansForStatement -= 1;
                }
              }
              nonZerosLeft -= jacobiansForStatement;  // Update nonzeros so that we know if it is the last round.

              Identifier storedIdentifier = lastIdentifier;
              // storeManual creates a new identifier which is either the identifier of the output w or the temporary
              // staggering variables t_1, t_2, ...
              tape.storeManual(value.getValue(), lastIdentifier, jacobiansForStatement + (int)staggeringActive);
              if (staggeringActive) {  // Not the first staggering so push the last output.
                tape.pushJacobiManual(1.0, 0.0, storedIdentifier);
              }

              // Push the rest of the Jacobians for the statement.
              while (jacobiansForStatement > 0) {
                if (Real() != (Real)jacobiean(curOut, curIn)) {
                  tape.pushJacobiManual(jacobiean(curOut, curIn), 0.0, inputData[curIn]);
                  jacobiansForStatement -= 1;
                }
                curIn += 1;
              }

              staggeringActive = true;
            }

            value.getIdentifier() = lastIdentifier; /* now set gradient data for the real output value */
          } else {
            // Disable tape index since there is no dependency.
            tape.destroyIdentifier(value.value(), value.getIdentifier());
          }
        }
      }
  };

#ifndef DOXYGEN_DISABLE
  /**
   * @brief Helper implementation of the PreaccumulationHelper interface for forward AD tapes.
   *
   * This implementation does nothing in all methods.
   */
  struct PreaccumulationHelperNoOpBase {
    public:

      /// Does nothing.
      template<typename... Inputs>
      void addInput(Inputs const&... inputs) {
        CODI_UNUSED(inputs...);
        // Do nothing.
      }

      /// Does nothing.
      template<typename... Inputs>
      void start(Inputs const&... inputs) {
        CODI_UNUSED(inputs...);
        // Do nothing.
      }

      /// Does nothing.
      template<typename... Outputs>
      void addOutput(Outputs&... outputs) {
        CODI_UNUSED(outputs...);
        // Do nothing.
      }

      /// Does nothing.
      template<typename... Outputs>
      void finish(bool const storeAdjoints, Outputs&... outputs) {
        CODI_UNUSED(storeAdjoints, outputs...);
        // Do nothing.
      }
  };

  /// Specialize PreaccumulationHelper for forward tapes.
  template<typename Type>
  struct PreaccumulationHelper<Type, TapeTraits::EnableIfForwardTape<typename Type::Tape>>
      : public PreaccumulationHelperNoOpBase {};

  /// Specialize PreaccumulationHelper for doubles.
  template<>
  struct PreaccumulationHelper<double, void> : public PreaccumulationHelperNoOpBase {};
#endif
}
