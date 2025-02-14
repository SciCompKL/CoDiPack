/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2025 Chair for Scientific Computing (SciComp), University of Kaiserslautern-Landau
 * Homepage: http://scicomp.rptu.de
 * Contact:  Prof. Nicolas R. Gauger (codi@scicomp.uni-kl.de)
 *
 * Lead developers: Max Sagebaum, Johannes Blühdorn (SciComp, University of Kaiserslautern-Landau)
 *
 * This file is part of CoDiPack (http://scicomp.rptu.de/software/codi).
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

#include <vector>

#include "../../config.h"
#include "../../expressions/lhsExpressionInterface.hpp"
#include "../../misc/exceptions.hpp"
#include "../../tapes/interfaces/fullTapeInterface.hpp"
#include "../../tapes/tagging/tagTapeReverse.hpp"
#include "../../traits/gradientTraits.hpp"
#include "../../traits/tapeTraits.hpp"
#include "../algorithms.hpp"
#include "../data/customAdjoints.hpp"
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
      using Type = CODI_DD(T_Type, CODI_DEFAULT_LHS_EXPRESSION);

      using Real = typename Type::Real;              ///< See LhsExpressionInterface.
      using Identifier = typename Type::Identifier;  ///< See LhsExpressionInterface.
      using Gradient = typename Type::Gradient;      ///< See LhsExpressionInterface.

      /// See LhsExpressionInterface.
      using Tape = CODI_DD(typename Type::Tape, CODI_DEFAULT_TAPE);
      using Position = typename Tape::Position;  ///< See PositionalEvaluationTapeInterface.

      std::vector<Identifier> inputData;   ///< List of input identifiers. Can be added manually after start() was
                                           ///< called.
      std::vector<Identifier> outputData;  ///< List of output identifiers. Can be added manually before finish() is
                                           ///< called. Has to be in sync with outputValues.
      std::vector<Type*> outputValues;     ///< List of output value pointers. Can be added manually before finish() is
                                           ///< called. Has to be in sync with outputData.

      std::vector<Gradient> localAdjoints;  ///< Vector of local adjoint variables. Persists across preaccumulations to
                                            ///< reduce the number of allocations, can be freed anytime if needed.

    protected:

      Position startPos;                        ///< Starting position for the region.
      std::vector<Gradient> storedAdjoints;     ///< If adjoints of inputs should be stored, before the preaccumulation.
      JacobianCountNonZerosRow<Real> jacobian;  ///< Jacobian for the preaccumulation.

    public:

      /// Constructor
      PreaccumulationHelper()
          : inputData(), outputData(), outputValues(), localAdjoints(), startPos(), storedAdjoints(), jacobian(0, 0) {}

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

        EventSystem<Tape>::notifyPreaccStartListeners(tape);

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

    protected:
      /// Internal implementation of workflow for preaccumulation with local adjoints.
      template<typename Func, typename... Outputs>
      void finishInternal(Func& coreRoutine, Outputs&... outputs) {
        Tape& tape = Type::getTape();

        if (tape.isActive()) {
          addOutputRecursive(outputs...);

          tape.setPassive();
          coreRoutine();
          storeJacobian();
          tape.setActive();
        }

        EventSystem<Tape>::notifyPreaccFinishListeners(tape);
      }

    public:
      /// Finish the preaccumulation region and perform the preaccumulation. See `addOutput()` for outputs.
      /// Not compatible with simultaneous thread-local preaccumulations with shared inputs. In this case, see
      /// finishLocalMappedAdjoints, finishLocalAdjointVectorPreprocessTape, finishLocalAdjoints,
      /// finishLocalAdjointVector, and finishLocalAdjointVectorOffset.
      template<typename... Outputs>
      void finish(bool const storeAdjoints, Outputs&... outputs) {
        auto coreRoutine = [&storeAdjoints, this]() {
          if (storeAdjoints) {
            storeInputAdjoints();
          }
          computeJacobian();
          if (storeAdjoints) {
            restoreInputAdjoints();
          }
        };

        finishInternal(coreRoutine, outputs...);
      }

      /// Finish the preaccumulation region and perform the preaccumulation. Creates a local map of adjoints instead of
      /// using adjoints from the tape. See `addOutput()` for outputs.
      template<typename... Outputs>
      void finishLocalMappedAdjoints(Outputs&... outputs) {
        auto coreRoutine = [this]() {
          computeJacobianLocalMappedAdjoints();
        };

        finishInternal(coreRoutine, outputs...);
      }

      /// Finish the preaccumulation region and perform the preaccumulation. Create a local adjoint vector instead of
      /// using adjoints from the tape. Edits the tape, remapping the identifiers to a contiguous range.
      /// More efficient than finishLocalMappedAdjoints if both the numbers of inputs and outputs are > 1. Behaves like
      /// finishLocalMappedAdjoints if the underlying tape does not support editing.
      template<typename... Outputs>
      void finishLocalAdjointVectorPreprocessTape(Outputs&... outputs) {
        auto coreRoutine = [this]() {
          computeJacobianLocalAdjointVectorPreprocessTapeIfAvailable<Tape>();  // otherwise
                                                                               // computeJacobianLocalMappedAdjoints
        };

        finishInternal(coreRoutine, outputs...);
      }

      /// Finish the preaccumulation region and perform the preaccumulation. Uses local adjoints instead of adjoints
      /// from the tape. Behaves either like finishLocalMappedAdjoints or like finishLocalAdjointVectorPreprocessTape,
      /// depending on which is more efficient given the numbers of inputs and outputs. See `addOutput()` for outputs.
      template<typename... Outputs>
      void finishLocalAdjoints(Outputs&... outputs) {
        auto coreRoutine = [this]() {
          if (std::min(inputData.size(), outputData.size()) > 1) {
            computeJacobianLocalAdjointVectorPreprocessTapeIfAvailable<Tape>();  // otherwise
                                                                                 // computeJacobianLocalMappedAdjoints
          } else {
            computeJacobianLocalMappedAdjoints();
          }
        };

        finishInternal(coreRoutine, outputs...);
      }

      /// Finish the preaccumulation region and perform the preaccumulation. Maintains a local adjoint vector that is as
      /// large as the global one. See `addOutput()` for outputs.
      template<typename... Outputs>
      void finishLocalAdjointVector(Outputs&... outputs) {
        auto coreRoutine = [this]() {
          computeJacobianLocalAdjointVector();
        };

        finishInternal(coreRoutine, outputs...);
      }

      /// Finish the preaccumulation region and perform the preaccumulation. Precomputes the identifier range used in
      /// the recording and creates a local adjoint vector, into which we address with an offset. Depends on tape
      /// editing features to preprocess identifiers. If tape editing is not supported, we fall back to the behaviour of
      /// finishLocalAdjointVector(). See `addOutput()` for outputs.
      template<typename... Outputs>
      void finishLocalAdjointVectorOffset(Outputs&... outputs) {
        auto coreRoutine = [this]() {
          computeJacobianLocalAdjointVectorOffsetIfAvailable<Tape>();  // otherwise
                                                                       // computeJacobianLocalAdjointVector
        };

        finishInternal(coreRoutine, outputs...);
      }

    private:

      // Tape supports editing -> use a map to edit its identifiers. Disabled by SFINAE otherwise.
      template<typename Tape>
      TapeTraits::EnableIfSupportsEditing<Tape> computeJacobianLocalAdjointVectorPreprocessTapeIfAvailable() {
        computeJacobianLocalAdjointVectorPreprocessTape();
      }

      // Tape does not support editing -> use a map for the adjoints. Disabled by SFINAE otherwise.
      template<typename Tape>
      TapeTraits::EnableIfNoEditing<Tape> computeJacobianLocalAdjointVectorPreprocessTapeIfAvailable() {
        computeJacobianLocalMappedAdjoints();
      }

      // Tape supports editing -> use a map to edit its identifiers. Disabled by SFINAE otherwise.
      template<typename Tape>
      TapeTraits::EnableIfSupportsEditing<Tape> computeJacobianLocalAdjointVectorOffsetIfAvailable() {
        computeJacobianLocalAdjointVectorOffset();
      }

      // Tape does not support editing -> use a map for the adjoints. Disabled by SFINAE otherwise.
      template<typename Tape>
      TapeTraits::EnableIfNoEditing<Tape> computeJacobianLocalAdjointVectorOffsetIfAvailable() {
        computeJacobianLocalAdjointVector();
      }

      void addInputLogic(Type const& input) {
        EventSystem<Tape>::notifyPreaccAddInputListeners(Type::getTape(), input.getValue(), input.getIdentifier());
        Identifier const& identifier = input.getIdentifier();
        if (Type::getTape().getPassiveIndex() != identifier) {
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
        EventSystem<Tape>::notifyPreaccAddOutputListeners(Type::getTape(), output.value(), output.getIdentifier());
        Identifier const& identifier = output.getIdentifier();
        if (Type::getTape().getPassiveIndex() != identifier) {
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

      void resizeJacobian() {
        if (jacobian.getM() != outputData.size() || jacobian.getN() != inputData.size()) {
          jacobian.resize(outputData.size(), inputData.size());
        }
      }

      void computeJacobian() {
        // Perform the accumulation of the tape part.
        Tape& tape = Type::getTape();
        Position endPos = tape.getPosition();

        resizeJacobian();

        // Manage adjoints manually to reduce the impact of locking on the performance.
        tape.resizeAdjointVector();
        tape.beginUseAdjointVector();

        Algorithms<Type, false>::computeJacobian(startPos, endPos, inputData.data(), inputData.size(),
                                                 outputData.data(), outputData.size(), jacobian,
                                                 AdjointsManagement::Manual);

        tape.resetTo(startPos, true, AdjointsManagement::Manual);

        tape.endUseAdjointVector();
      }

      void computeJacobianLocalAdjointVector() {
        // Perform the accumulation of the tape part.
        Tape& tape = Type::getTape();
        Position endPos = tape.getPosition();

        resizeJacobian();

        size_t requiredVectorSize = tape.getParameter(TapeParameters::LargestIdentifier) + 1;

        this->localAdjoints.resize(requiredVectorSize);

        Algorithms<Type, false>::computeJacobianCustomAdjoints(startPos, endPos, inputData.data(), inputData.size(),
                                                               outputData.data(), outputData.size(), jacobian,
                                                               this->localAdjoints.data());

        tape.resetTo(startPos, false);
      }

      void computeJacobianLocalAdjointVectorOffset() {
        // Perform the accumulation of the tape part.
        Tape& tape = Type::getTape();
        Position endPos = tape.getPosition();

        resizeJacobian();

        // Determine minimum and maximum identifier used in the recording.

        Identifier minIdentifier = std::numeric_limits<Identifier>::max();
        Identifier maxIdentifier = std::numeric_limits<Identifier>::min();

        auto determineMinMaxIdentifier = [&minIdentifier, &maxIdentifier](typename Tape::Identifier const& identifier) {
          minIdentifier = std::min(minIdentifier, identifier);
          maxIdentifier = std::max(maxIdentifier, identifier);
        };

        // Begin by processing inputs and outputs.
        for (auto const& identifier : inputData) {
          determineMinMaxIdentifier(identifier);
        }

        for (auto const& identifier : outputData) {
          determineMinMaxIdentifier(identifier);
        }

        // Process the tape. Does not edit identifiers in the tape.
        tape.template editIdentifiers(determineMinMaxIdentifier, startPos, endPos);

        // Plus one to cover the range [minIdentifier, maxIdentifier].
        size_t requiredVectorSize = maxIdentifier - minIdentifier + 1;
        this->localAdjoints.resize(requiredVectorSize);

        // Define adjoints that take into account the offset when addressing into the vector.
        using LocalAdjointsOffset = AdjointVectorWithOffset<Identifier, Gradient>;
        LocalAdjointsOffset localAdjointsOffset(this->localAdjoints.data(), minIdentifier);

        // Preaccumulation with a local adjoint vector and identifier offsets.
        Algorithms<Type, false>::computeJacobianCustomAdjoints(startPos, endPos, inputData.data(), inputData.size(),
                                                               outputData.data(), outputData.size(), jacobian,
                                                               localAdjointsOffset);

        tape.resetTo(startPos, false);
      }

      void computeJacobianLocalMappedAdjoints() {
        // Perform the accumulation of the tape part.
        Tape& tape = Type::getTape();
        Position endPos = tape.getPosition();

        resizeJacobian();

        // Create a local map with adjoints.
        using LocalMappedAdjoints = MappedAdjoints<typename Tape::Identifier, typename Tape::Gradient>;
        LocalMappedAdjoints mappedAdjoints;

        Algorithms<Type, false>::computeJacobianCustomAdjoints(startPos, endPos, inputData.data(), inputData.size(),
                                                               outputData.data(), outputData.size(), jacobian,
                                                               mappedAdjoints);

        tape.resetTo(startPos, false);
      }

      void computeJacobianLocalAdjointVectorPreprocessTape() {
        // Perform the accumulation of the tape part.
        Tape& tape = Type::getTape();
        Position endPos = tape.getPosition();

        resizeJacobian();

        // Used internally for remapping identifiers.
        using IdentifierMap = std::map<typename Tape::Identifier, typename Tape::Identifier>;

        // Build a map of identifiers, remapping identifiers in the recording to contiguous ones.
        auto nextIdentifier = typename Tape::Identifier() + 1;
        IdentifierMap oldToNewIdentifierMap;

        // If needed, inserts the old identifier into the map and associates it with the next identifier. Either way,
        // returns the associated new identifier.
        auto accessOldToNewIdentifierMap = [&](typename Tape::Identifier const& oldIdentifier) -> typename Tape::Identifier const& {
          auto result = oldToNewIdentifierMap.insert({oldIdentifier, nextIdentifier});
          if (result.second) {  // insertion took place
            ++nextIdentifier;
          }
          return result.first->second;
        };

        // Remap input identifiers explicitly to account for inputs that are actually not used in the recording.
        for (auto const& oldIdentifier : inputData) {
          accessOldToNewIdentifierMap(oldIdentifier);
        }

        // Remap output identifiers explicitly to account for outputs that actually do not depend on the inputs.
        for (auto const& oldIdentifier : outputData) {
          accessOldToNewIdentifierMap(oldIdentifier);
        }

        auto editIdentifier = [&](typename Tape::Identifier& oldIdentifier) {
          oldIdentifier = accessOldToNewIdentifierMap(oldIdentifier);
        };

        // Process the recording to complete the map, edit the tape on the fly.
        tape.template editIdentifiers(editIdentifier, startPos, endPos);

        // Build new vectors of input and output identifiers.
        std::vector<typename Tape::Identifier> newInputData;
        newInputData.reserve(inputData.size());
        for (auto const& identifier : inputData) {
          newInputData.push_back(oldToNewIdentifierMap[identifier]);
        }

        std::vector<typename Tape::Identifier> newOutputData;
        newOutputData.reserve(outputData.size());
        for (auto const& identifier : outputData) {
          newOutputData.push_back(oldToNewIdentifierMap[identifier]);
        }

        // The association with the original input/output identifiers and Jacobian entries is made by position, so we no
        // longer need the identifier map.
        oldToNewIdentifierMap.clear();

        // Create local adjoints. nextIdentifier holds the local adjoint vector size.
        std::vector<typename Tape::Gradient> localAdjoints(nextIdentifier);

        // Preaccumulation with remapped identifiers on local adjoints.
        Algorithms<Type, false>::computeJacobianCustomAdjoints(startPos, endPos, newInputData.data(),
                                                               newInputData.size(), newOutputData.data(),
                                                               newOutputData.size(), jacobian, localAdjoints.data());

        tape.resetTo(startPos, false);
      }

      void storeJacobian() {
        Tape& tape = Type::getTape();

        for (size_t curOut = 0; curOut < outputData.size(); ++curOut) {
          Type& value = *outputValues[curOut];
          if (0 != jacobian.nonZerosRow(curOut)) {
            int nonZerosLeft = jacobian.nonZerosRow(curOut);
            jacobian.nonZerosRow(curOut) = 0;

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
                tape.pushJacobianManual(1.0, 0.0, storedIdentifier);
              }

              // Push the rest of the Jacobians for the statement.
              while (jacobiansForStatement > 0) {
                if (Real() != (Real)jacobian(curOut, curIn)) {
                  tape.pushJacobianManual(jacobian(curOut, curIn), 0.0, inputData[curIn]);
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

      /// Does nothing.
      template<typename... Outputs>
      void finishLocalMappedAdjoints(Outputs&... outputs) {
        CODI_UNUSED(outputs...);
        // Do nothing.
      }

      /// Does nothing.
      template<typename... Outputs>
      void finishLocalAdjointVectorPreprocessTape(Outputs&... outputs) {
        CODI_UNUSED(outputs...);
        // Do nothing.
      }

      /// Does nothing.
      template<typename... Outputs>
      void finishLocalAdjoints(Outputs&... outputs) {
        CODI_UNUSED(outputs...);
        // Do nothing.
      }

      /// Does nothing.
      template<typename... Outputs>
      void finishLocalAdjointVector(Outputs&... outputs) {
        CODI_UNUSED(outputs...);
        // Do nothing.
      }

      /// Does nothing.
      template<typename... Outputs>
      void finishLocalAdjointVectorOffset(Outputs&... outputs) {
        CODI_UNUSED(outputs...);
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

  /**
   * @brief Helper implementation of the PreaccumulationHelper interface for tag tapes.
   *
   * It changes the current tag of the tape and sets this tag on all inputs. After the
   * preaccumulation is finished the original tag is restored on the tape, the inputs
   * and the outputs.
   */
  template<typename T_Type>
  struct PreaccumulationHelper<T_Type, TapeTraits::EnableIfTagTapeReverse<typename T_Type::Tape>> {
    public:

      using Type = CODI_DD(T_Type, CODI_DEFAULT_LHS_EXPRESSION);  ///< See PreaccumulationHelper.

      using Tape = CODI_DD(typename Type::Tape, CODI_T(TagTapeReverse<Type, int>));  ///< Tape of the CoDiPack type.
      using Tag = typename Tape::Tag;                                                ///< Tag of the tape.

    private:

      std::vector<Type const*> inputLocations;
      std::vector<Type*> outputLocations;
      Tag oldTag;

    public:

      /// Constructor.
      CODI_INLINE PreaccumulationHelper() = default;

      /// Gathers the input values.
      template<typename... Inputs>
      void addInput(Inputs const&... inputs) {
        Tape& tape = getTape();

        if (tape.isActive() && tape.isPreaccumulationHandlingEnabled()) {
          addInputRecursive(inputs...);
        }
      }

      /// Set special tag on inputs.
      template<typename... Inputs>
      void start(Inputs const&... inputs) {
        Tape& tape = getTape();

        if (tape.isActive() && tape.isPreaccumulationHandlingEnabled()) {
          inputLocations.clear();
          outputLocations.clear();
          oldTag = tape.getCurTag();
          tape.setCurTag(tape.getPreaccumulationHandlingTag());

          addInputRecursive(inputs...);
        }
      }

      /// Gathers the outputs.
      template<typename... Outputs>
      void addOutput(Outputs&... outputs) {
        Tape& tape = getTape();

        if (tape.isActive() && tape.isPreaccumulationHandlingEnabled()) {
          addOutputRecursive(outputs...);
        }
      }

      /// Reverts the tags on all input and output values.
      template<typename... Outputs>
      void finish(bool const storeAdjoints, Outputs&... outputs) {
        CODI_UNUSED(storeAdjoints);

        Tape& tape = getTape();

        if (tape.isActive() && tape.isPreaccumulationHandlingEnabled()) {
          addOutputRecursive(outputs...);

          tape.setCurTag(oldTag);
          for (Type const* curInput : inputLocations) {
            tape.setTagOnVariable(*curInput);
          }
          for (Type* curOutput : outputLocations) {
            tape.setTagOnVariable(*curOutput);
          }
        }
      }

      /// Reverts the tags on all input and output values.
      template<typename... Outputs>
      void finishLocalMappedAdjoints(Outputs&... outputs) {
        finish(false, outputs...);
      }

      /// Reverts the tags on all input and output values.
      template<typename... Outputs>
      void finishLocalAdjointVectorPreprocessTape(Outputs&... outputs) {
        finish(false, outputs...);
      }

      /// Reverts the tags on all input and output values.
      template<typename... Outputs>
      void finishLocalAdjoints(Outputs&... outputs) {
        finish(false, outputs...);
      }

      /// Reverts the tags on all input and output values.
      template<typename... Outputs>
      void finishLocalAdjointVector(Outputs&... outputs) {
        finish(false, outputs...);
      }

      /// Reverts the tags on all input and output values.
      template<typename... Outputs>
      void finishLocalAdjointVectorOffset(Outputs&... outputs) {
        finish(false, outputs...);
      }

    private:

      /// Terminator for the recursive implementation.
      void addInputRecursive() {
        // Terminator implementation.
      }

      template<typename... Inputs>
      void addInputRecursive(Type const& input, Inputs const&... r) {
        handleInput(input);
        addInputRecursive(r...);
      }

      void handleInput(Type const& input) {
        if (Type::getTape().getPassiveIndex() != input.getIdentifier()) {
          inputLocations.push_back(&input);
          getTape().setTagOnVariable(input);
        }
      }

      /// Terminator for the recursive implementation.
      void addOutputRecursive() {
        // Terminator implementation.
      }

      template<typename... Outputs>
      void addOutputRecursive(Type& output, Outputs&... r) {
        handleOutput(output);
        addOutputRecursive(r...);
      }

      void handleOutput(Type& value) {
        outputLocations.push_back(&value);
      }

      Tape& getTape() {
        return Type::getTape();
      }
  };
#endif
}
