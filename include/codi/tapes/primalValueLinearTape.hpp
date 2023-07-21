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

#include <algorithm>
#include <functional>
#include <type_traits>

#include "../config.h"
#include "../expressions/lhsExpressionInterface.hpp"
#include "../expressions/logic/compileTimeTraversalLogic.hpp"
#include "../expressions/logic/constructStaticContext.hpp"
#include "../expressions/logic/traversalLogic.hpp"
#include "../misc/macros.hpp"
#include "../misc/memberStore.hpp"
#include "../traits/expressionTraits.hpp"
#include "data/chunk.hpp"
#include "indices/indexManagerInterface.hpp"
#include "primalValueBaseTape.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Final implementation for a primal value tape with a linear index management.
   *
   * This class implements the interface methods from the PrimalValueBaseTape.
   *
   * @tparam T_TapeTypes  JacobianTapeTypes definition.
   */
  template<typename T_TapeTypes>
  struct PrimalValueLinearTape : public PrimalValueBaseTape<T_TapeTypes, PrimalValueLinearTape<T_TapeTypes>> {
    public:

      using TapeTypes =
          CODI_DD(T_TapeTypes,
                  CODI_T(PrimalValueTapeTypes<double, double, IndexManagerInterface<int>, StatementEvaluatorInterface,
                                              DefaultChunkedData>));  ///< See PrimalValueLinearTape.

      using Base = PrimalValueBaseTape<T_TapeTypes, PrimalValueLinearTape<T_TapeTypes>>;  ///< Base class abbreviation.
      friend Base;  ///< Allow the base class to call protected and private methods.

      using Real = typename TapeTypes::Real;                              ///< See TapeTypesInterface.
      using Gradient = typename TapeTypes::Gradient;                      ///< See TapeTypesInterface.
      using IndexManager = typename TapeTypes::IndexManager;              ///< See PrimalValueTapeTypes.
      using Identifier = typename TapeTypes::Identifier;                  ///< See TapeTypesInterface.
      using PassiveReal = RealTraits::PassiveReal<Real>;                  ///< Basic computation type.
      using StatementEvaluator = typename TapeTypes::StatementEvaluator;  ///< See PrimalValueTapeTypes.
      using EvalHandle = typename TapeTypes::EvalHandle;                  ///< See PrimalValueTapeTypes.
      using Position = typename Base::Position;                           ///< See TapeTypesInterface.

      /// Constructor
      PrimalValueLinearTape() : Base() {}

      using Base::clearAdjoints;

      /// \copydoc codi::PositionalEvaluationTapeInterface::clearAdjoints
      /// <br> Implementation: Automatic adjoints management has no effect. Primal value tapes do not implement adjoints
      /// locking.
      void clearAdjoints(Position const& start, Position const& end,
                         AdjointsManagement adjointsManagement = AdjointsManagement::Automatic) {
        CODI_UNUSED(adjointsManagement);

        using IndexPosition = CODI_DD(typename IndexManager::Position, int);
        IndexPosition startIndex = this->externalFunctionData.template extractPosition<IndexPosition>(start);
        IndexPosition endIndex = this->externalFunctionData.template extractPosition<IndexPosition>(end);

        startIndex = std::min(startIndex, (IndexPosition)this->adjoints.size() - 1);
        endIndex = std::min(endIndex, (IndexPosition)this->adjoints.size() - 1);

        for (IndexPosition curPos = endIndex + 1; curPos <= startIndex; curPos += 1) {
          this->adjoints[curPos] = Gradient();
        }
      }

    protected:

      /// \copydoc codi::PrimalValueBaseTape::internalEvaluateForward_Step3_EvalStatements
      CODI_INLINE static void internalEvaluateForward_Step3_EvalStatements(
          /* data from call */
          PrimalValueLinearTape& tape, Real* primalVector, ADJOINT_VECTOR_TYPE* adjointVector,
          /* data from constantValueData */
          size_t& curConstantPos, size_t const& endConstantPos, PassiveReal const* const constantValues,
          /* data from passiveValueData */
          size_t& curPassivePos, size_t const& endPassivePos, Real const* const passiveValues,
          /* data from rhsIdentifiersData */
          size_t& curRhsIdentifiersPos, size_t const& endRhsIdentifiersPos, Identifier const* const rhsIdentifiers,
          /* data from statementData */
          size_t& curStatementPos, size_t const& endStatementPos,
          Config::ArgumentSize const* const numberOfPassiveArguments, EvalHandle const* const stmtEvalhandle,
          /* data from index handler */
          size_t const& startAdjointPos, size_t const& endAdjointPos) {
        CODI_UNUSED(endConstantPos, endPassivePos, endRhsIdentifiersPos, endStatementPos);

        size_t curAdjointPos = startAdjointPos;

        while (curAdjointPos < endAdjointPos) CODI_Likely {
          curAdjointPos += 1;

          Config::ArgumentSize nPassiveValues = numberOfPassiveArguments[curStatementPos];

          if (Config::StatementInputTag != nPassiveValues) CODI_Likely {
            Gradient lhsTangent = Gradient();

            primalVector[curAdjointPos] = StatementEvaluator::template callForward<PrimalValueLinearTape>(
                stmtEvalhandle[curStatementPos], primalVector, adjointVector, lhsTangent, nPassiveValues,
                curConstantPos, constantValues, curPassivePos, passiveValues, curRhsIdentifiersPos, rhsIdentifiers);

#if CODI_VariableAdjointInterfaceInPrimalTapes
            adjointVector->setLhsTangent(curAdjointPos);
            EventSystem<PrimalValueLinearTape>::notifyStatementEvaluateListeners(
                tape, curAdjointPos, adjointVector->getVectorSize(), adjointVector->getAdjointVec(curAdjointPos));
#else
            adjointVector[curAdjointPos] = lhsTangent;

            EventSystem<PrimalValueLinearTape>::notifyStatementEvaluateListeners(
                tape, curAdjointPos, GradientTraits::dim<Gradient>(), GradientTraits::toArray(lhsTangent).data());
#endif
            EventSystem<PrimalValueLinearTape>::notifyStatementEvaluatePrimalListeners(tape, curAdjointPos,
                                                                                       primalVector[curAdjointPos]);
          }

          curStatementPos += 1;
        }
      }

      /// \copydoc codi::PrimalValueBaseTape::internalEvaluatePrimal_Step3_EvalStatements
      CODI_INLINE static void internalEvaluatePrimal_Step3_EvalStatements(
          /* data from call */
          PrimalValueLinearTape& tape, Real* primalVector,
          /* data from constantValueData */
          size_t& curConstantPos, size_t const& endConstantPos, PassiveReal const* const constantValues,
          /* data from passiveValueData */
          size_t& curPassivePos, size_t const& endPassivePos, Real const* const passiveValues,
          /* data from rhsIdentifiersData */
          size_t& curRhsIdentifiersPos, size_t const& endRhsIdentifiersPos, Identifier const* const rhsIdentifiers,
          /* data from statementData */
          size_t& curStatementPos, size_t const& endStatementPos,
          Config::ArgumentSize const* const numberOfPassiveArguments, EvalHandle const* const stmtEvalhandle,
          /* data from index handler */
          size_t const& startAdjointPos, size_t const& endAdjointPos) {
        CODI_UNUSED(endConstantPos, endPassivePos, endRhsIdentifiersPos, endStatementPos);

        size_t curAdjointPos = startAdjointPos;

        while (curAdjointPos < endAdjointPos) CODI_Likely {
          curAdjointPos += 1;

          Config::ArgumentSize nPassiveValues = numberOfPassiveArguments[curStatementPos];

          if (Config::StatementInputTag != nPassiveValues) CODI_Likely {
            primalVector[curAdjointPos] = StatementEvaluator::template callPrimal<PrimalValueLinearTape>(
                stmtEvalhandle[curStatementPos], primalVector, nPassiveValues, curConstantPos, constantValues,
                curPassivePos, passiveValues, curRhsIdentifiersPos, rhsIdentifiers);

            EventSystem<PrimalValueLinearTape>::notifyStatementEvaluatePrimalListeners(tape, curAdjointPos,
                                                                                       primalVector[curAdjointPos]);
          }

          curStatementPos += 1;
        }
      }

      /// \copydoc codi::PrimalValueBaseTape::internalEvaluateReverse_Step3_EvalStatements
      CODI_INLINE static void internalEvaluateReverse_Step3_EvalStatements(
          /* data from call */
          PrimalValueLinearTape& tape, Real* primalVector, ADJOINT_VECTOR_TYPE* adjointVector,
          /* data from constantValueData */
          size_t& curConstantPos, size_t const& endConstantPos, PassiveReal const* const constantValues,
          /* data from passiveValueData */
          size_t& curPassivePos, size_t const& endPassivePos, Real const* const passiveValues,
          /* data from rhsIdentifiersData */
          size_t& curRhsIdentifiersPos, size_t const& endRhsIdentifiersPos, Identifier const* const rhsIdentifiers,
          /* data from statementData */
          size_t& curStatementPos, size_t const& endStatementPos,
          Config::ArgumentSize const* const numberOfPassiveArguments, EvalHandle const* const stmtEvalhandle,
          /* data from index handler */
          size_t const& startAdjointPos, size_t const& endAdjointPos) {
        CODI_UNUSED(endConstantPos, endPassivePos, endRhsIdentifiersPos, endStatementPos);

        size_t curAdjointPos = startAdjointPos;

        while (curAdjointPos > endAdjointPos) CODI_Likely {
          curStatementPos -= 1;

          Config::ArgumentSize nPassiveValues = numberOfPassiveArguments[curStatementPos];

          if (Config::StatementInputTag != nPassiveValues) CODI_Likely {
#if CODI_VariableAdjointInterfaceInPrimalTapes

            EventSystem<PrimalValueLinearTape>::notifyStatementEvaluateListeners(
                tape, curAdjointPos, adjointVector->getVectorSize(), adjointVector->getAdjointVec(curAdjointPos));

            Gradient const lhsAdjoint{};
            adjointVector->setLhsAdjoint(curAdjointPos);
#else
            Gradient const lhsAdjoint = adjointVector[curAdjointPos];

            EventSystem<PrimalValueLinearTape>::notifyStatementEvaluateListeners(
                tape, curAdjointPos, GradientTraits::dim<Gradient>(), GradientTraits::toArray(lhsAdjoint).data());

            if (Config::ReversalZeroesAdjoints) {
              adjointVector[curAdjointPos] = Gradient();
            }
#endif
            EventSystem<PrimalValueLinearTape>::notifyStatementEvaluatePrimalListeners(tape, curAdjointPos,
                                                                                       primalVector[curAdjointPos]);

            StatementEvaluator::template callReverse<PrimalValueLinearTape>(
                stmtEvalhandle[curStatementPos], primalVector, adjointVector, lhsAdjoint, nPassiveValues,
                curConstantPos, constantValues, curPassivePos, passiveValues, curRhsIdentifiersPos, rhsIdentifiers);
          }

          curAdjointPos -= 1;
        }
      }

      /// \copydoc codi::PrimalValueBaseTape::internalResetPrimalValues
      /// Empty implementation; primal values are not overwritten with linear index management.
      CODI_INLINE void internalResetPrimalValues(Position const& pos) {
        CODI_UNUSED(pos);

        // Nothing to do.
      }

      /// \copydoc codi::PrimalValueBaseTape::pushStmtData
      /// Only the number of arguments is required for linear index managers.
      CODI_INLINE void pushStmtData(Identifier const& index, Config::ArgumentSize const& numberOfPassiveArguments,
                                    Real const& oldPrimalValue, EvalHandle evalHandle) {
        CODI_UNUSED(index, oldPrimalValue);

        Base::statementData.pushData(numberOfPassiveArguments, evalHandle);
      }

    public:
      /// \copydoc codi::PrimalValueBaseTape::revertPrimals
      /// Empty implementation; primal values are not overwritten with linear index management.
      void revertPrimals(Position const& pos) {
        CODI_UNUSED(pos);

        // Primal values do not need to be reset.
      }
  };
}
