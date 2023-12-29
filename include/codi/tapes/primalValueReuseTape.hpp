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
   * @brief Final implementation for a primal value tape with a reuse index management.
   *
   * This class implements the interface methods from the PrimalValueBaseTape.
   *
   * @tparam T_TapeTypes  JacobianTapeTypes definition.
   */
  template<typename T_TapeTypes>
  struct PrimalValueReuseTape : public PrimalValueBaseTape<T_TapeTypes, PrimalValueReuseTape<T_TapeTypes>> {
    public:

      using TapeTypes =
          CODI_DD(T_TapeTypes,
                  CODI_T(PrimalValueTapeTypes<double, double, IndexManagerInterface<int>, StatementEvaluatorInterface,
                                              DefaultChunkedData>));  ///< See PrimalValueReuseTape.

      using Base = PrimalValueBaseTape<TapeTypes, PrimalValueReuseTape<TapeTypes>>;  ///< Base class abbreviation.
      friend Base;  ///< Allow the base class to call protected and private methods.

      using Real = typename TapeTypes::Real;                              ///< See TapeTypesInterface.
      using Gradient = typename TapeTypes::Gradient;                      ///< See TapeTypesInterface.
      using Identifier = typename TapeTypes::Identifier;                  ///< See TapeTypesInterface.
      using PassiveReal = RealTraits::PassiveReal<Real>;                  ///< Basic computation type.
      using StatementEvaluator = typename TapeTypes::StatementEvaluator;  ///< See PrimalValueTapeTypes.
      using EvalHandle = typename TapeTypes::EvalHandle;                  ///< See PrimalValueTapeTypes.
      using Position = typename Base::Position;                           ///< See TapeTypesInterface.

      using StatementData = typename TapeTypes::StatementData;  ///< See PrimalValueTapeTypes.

      /// Constructor
      PrimalValueReuseTape() : Base() {}

      using Base::clearAdjoints;

      /// \copydoc codi::PositionalEvaluationTapeInterface::clearAdjoints
      /// <br> Implementation: Automatic adjoints management has no effect. Primal value tapes do not implement adjoints
      /// locking.
      void clearAdjoints(Position const& start, Position const& end,
                         AdjointsManagement adjointsManagement = AdjointsManagement::Automatic) {
        CODI_UNUSED(adjointsManagement);

        auto clearFunc = [this](Identifier* lhsIndex, Config::ArgumentSize* passiveArgs, Real* oldPrimal,
                                EvalHandle* evalHandle) {
          CODI_UNUSED(passiveArgs, oldPrimal, evalHandle);

          if (*lhsIndex < (Identifier)this->adjoints.size()) {
            this->adjoints[*lhsIndex] = Gradient();
          }
        };

        using StmtPosition = typename StatementData::Position;
        StmtPosition startStmt = this->externalFunctionData.template extractPosition<StmtPosition>(start);
        StmtPosition endStmt = this->externalFunctionData.template extractPosition<StmtPosition>(end);

        this->statementData.forEachReverse(startStmt, endStmt, clearFunc);
      }

    protected:

      /// \copydoc codi::PrimalValueBaseTape::internalEvaluateForward_Step3_EvalStatements
      CODI_INLINE static void internalEvaluateForward_Step3_EvalStatements(
          /* data from call */
          PrimalValueReuseTape& tape, Real* primalVector, ADJOINT_VECTOR_TYPE* adjointVector,
          /* data from constantValueData */
          size_t& curConstantPos, size_t const& endConstantPos, PassiveReal const* const constantValues,
          /* data from passiveValueData */
          size_t& curPassivePos, size_t const& endPassivePos, Real const* const passiveValues,
          /* data from rhsIdentifiersData */
          size_t& curRhsIdentifiersPos, size_t const& endRhsIdentifiersPos, Identifier const* const rhsIdentifiers,
          /* data from statementData */
          size_t& curStatementPos, size_t const& endStatementPos, Identifier const* const lhsIdentifiers,
          Config::ArgumentSize const* const numberOfPassiveArguments, Real* const oldPrimalValues,
          EvalHandle const* const stmtEvalhandle) {
        CODI_UNUSED(endConstantPos, endPassivePos, endRhsIdentifiersPos);

        while (curStatementPos < endStatementPos) CODI_Likely {
          Identifier const lhsIdentifier = lhsIdentifiers[curStatementPos];

          Gradient lhsTangent = Gradient();

          oldPrimalValues[curStatementPos] = primalVector[lhsIdentifier];
          primalVector[lhsIdentifier] = StatementEvaluator::template callForward<PrimalValueReuseTape>(
              stmtEvalhandle[curStatementPos], primalVector, adjointVector, lhsTangent,
              numberOfPassiveArguments[curStatementPos], curConstantPos, constantValues, curPassivePos, passiveValues,
              curRhsIdentifiersPos, rhsIdentifiers);

#if CODI_VariableAdjointInterfaceInPrimalTapes
          adjointVector->setLhsTangent(lhsIdentifier);
          EventSystem<PrimalValueReuseTape>::notifyStatementEvaluateListeners(
              tape, lhsIdentifier, adjointVector->getVectorSize(), adjointVector->getAdjointVec(lhsIdentifier));
#else
          adjointVector[lhsIdentifier] = lhsTangent;
          EventSystem<PrimalValueReuseTape>::notifyStatementEvaluateListeners(
              tape, lhsIdentifier, GradientTraits::dim<Gradient>(), GradientTraits::toArray(lhsTangent).data());
#endif
          EventSystem<PrimalValueReuseTape>::notifyStatementEvaluatePrimalListeners(tape, lhsIdentifier,
                                                                                    primalVector[lhsIdentifier]);

          curStatementPos += 1;
        }
      }

      /// \copydoc codi::PrimalValueBaseTape::internalEvaluatePrimal_Step3_EvalStatements
      CODI_INLINE static void internalEvaluatePrimal_Step3_EvalStatements(
          /* data from call */
          PrimalValueReuseTape& tape, Real* primalVector,
          /* data from constantValueData */
          size_t& curConstantPos, size_t const& endConstantPos, PassiveReal const* const constantValues,
          /* data from passiveValueData */
          size_t& curPassivePos, size_t const& endPassivePos, Real const* const passiveValues,
          /* data from rhsIdentifiersData */
          size_t& curRhsIdentifiersPos, size_t const& endRhsIdentifiersPos, Identifier const* const rhsIdentifiers,
          /* data from statementData */
          size_t& curStatementPos, size_t const& endStatementPos, Identifier const* const lhsIdentifiers,
          Config::ArgumentSize const* const numberOfPassiveArguments, Real* const oldPrimalValues,
          EvalHandle const* const stmtEvalhandle) {
        CODI_UNUSED(endConstantPos, endPassivePos, endRhsIdentifiersPos);

        while (curStatementPos < endStatementPos) CODI_Likely {
          Identifier const lhsIdentifier = lhsIdentifiers[curStatementPos];

          oldPrimalValues[curStatementPos] = primalVector[lhsIdentifier];
          primalVector[lhsIdentifier] = StatementEvaluator::template callPrimal<PrimalValueReuseTape>(
              stmtEvalhandle[curStatementPos], primalVector, numberOfPassiveArguments[curStatementPos], curConstantPos,
              constantValues, curPassivePos, passiveValues, curRhsIdentifiersPos, rhsIdentifiers);

          EventSystem<PrimalValueReuseTape>::notifyStatementEvaluatePrimalListeners(tape, lhsIdentifier,
                                                                                    primalVector[lhsIdentifier]);

          curStatementPos += 1;
        }
      }

      /// \copydoc codi::PrimalValueBaseTape::internalEvaluateReverse_Step3_EvalStatements
      CODI_INLINE static void internalEvaluateReverse_Step3_EvalStatements(
          /* data from call */
          PrimalValueReuseTape& tape, Real* primalVector, ADJOINT_VECTOR_TYPE* adjointVector,
          /* data from constantValueData */
          size_t& curConstantPos, size_t const& endConstantPos, PassiveReal const* const constantValues,
          /* data from passiveValueData */
          size_t& curPassivePos, size_t const& endPassivePos, Real const* const passiveValues,
          /* data from rhsIdentifiersData */
          size_t& curRhsIdentifiersPos, size_t const& endRhsIdentifiersPos, Identifier const* const rhsIdentifiers,
          /* data from statementData */
          size_t& curStatementPos, size_t const& endStatementPos, Identifier const* const lhsIdentifiers,
          Config::ArgumentSize const* const numberOfPassiveArguments, Real const* const oldPrimalValues,
          EvalHandle const* const stmtEvalhandle) {
        CODI_UNUSED(endConstantPos, endPassivePos, endRhsIdentifiersPos);

        while (curStatementPos > endStatementPos) CODI_Likely {
          curStatementPos -= 1;

          Identifier const lhsIdentifier = lhsIdentifiers[curStatementPos];

#if CODI_VariableAdjointInterfaceInPrimalTapes
          EventSystem<PrimalValueReuseTape>::notifyStatementEvaluateListeners(
              tape, lhsIdentifier, adjointVector->getVectorSize(), adjointVector->getAdjointVec(lhsIdentifier));
          Gradient const lhsAdjoint{};
          adjointVector->setLhsAdjoint(lhsIdentifier);
#else
          Gradient const lhsAdjoint = adjointVector[lhsIdentifier];
          EventSystem<PrimalValueReuseTape>::notifyStatementEvaluateListeners(
              tape, lhsIdentifier, GradientTraits::dim<Gradient>(), GradientTraits::toArray(lhsAdjoint).data());
          adjointVector[lhsIdentifier] = Gradient();
#endif
          EventSystem<PrimalValueReuseTape>::notifyStatementEvaluatePrimalListeners(tape, lhsIdentifier,
                                                                                    primalVector[lhsIdentifier]);

          primalVector[lhsIdentifier] = oldPrimalValues[curStatementPos];

          StatementEvaluator::template callReverse<PrimalValueReuseTape>(
              stmtEvalhandle[curStatementPos], primalVector, adjointVector, lhsAdjoint,
              numberOfPassiveArguments[curStatementPos], curConstantPos, constantValues, curPassivePos, passiveValues,
              curRhsIdentifiersPos, rhsIdentifiers);
        }
      }

      /// \copydoc codi::PrimalValueBaseTape::internalResetPrimalValues
      CODI_INLINE void internalResetPrimalValues(Position const& pos) {
        // Reset primals.
        auto clearFunc = [this](Identifier* lhsIndex, Config::ArgumentSize* passiveArgs, Real* oldPrimal,
                                EvalHandle* evalHandle) {
          CODI_UNUSED(passiveArgs, evalHandle);

          this->primals[*lhsIndex] = *oldPrimal;
        };

        using StmtPosition = typename StatementData::Position;
        StmtPosition startStmt = this->externalFunctionData.template extractPosition<StmtPosition>(this->getPosition());
        StmtPosition endStmt = this->externalFunctionData.template extractPosition<StmtPosition>(pos);

        this->statementData.forEachReverse(startStmt, endStmt, clearFunc);
      }

      /// \copydoc codi::PrimalValueBaseTape::pushStmtData
      /// Only the number of arguments is required for linear index managers.
      CODI_INLINE void pushStmtData(Identifier const& index, Config::ArgumentSize const& numberOfPassiveArguments,
                                    Real const& oldPrimalValue, EvalHandle evalHandle) {
        Base::statementData.pushData(index, numberOfPassiveArguments, oldPrimalValue, evalHandle);
      }

    public:
      /// \copydoc codi::PrimalEvaluationTapeInterface::revertPrimals
      void revertPrimals(Position const& pos) {
        internalResetPrimalValues(pos);
      }
  };
}
