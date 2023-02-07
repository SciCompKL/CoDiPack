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

      using StmtFixedDataEntry = typename Base::StmtFixedDataEntry;  ///< See PrimalValueBaseTape.
      using StmtCallArgs = typename Base::StmtCallArgs;              ///< See PrimalValueBaseTape.

      template<typename T>
      using StackArray = typename Base::template StackArray<T>;  ///< See PrimalValueBaseTape.

      /// Constructor
      PrimalValueLinearTape() : Base() {}

      using Base::clearAdjoints;

      /// \copydoc codi::PositionalEvaluationTapeInterface::clearAdjoints
      void clearAdjoints(Position const& start, Position const& end) {
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
          /* data from dynamicSizeData */
          size_t& curDynamicSizePos, size_t const& endDynamicPos, char* const dynamicSizeValues,
          /* data from fixedSizeData */
          size_t& curFixedSizePos, size_t const& endFixedSizePos, char const* const fixedSizeValues,
          /* data from index handler */
          size_t const& startAdjointPos, size_t const& endAdjointPos) {
        CODI_UNUSED(endDynamicPos, endFixedSizePos);

        StackArray<Real> lhsPrimals{};
        StackArray<Gradient> lhsTangents{};
        StmtFixedDataEntry data;

        size_t curAdjointPos = startAdjointPos;

        while (curAdjointPos < endAdjointPos) {
          curFixedSizePos = data.readForward(fixedSizeValues, curFixedSizePos);

          if (Config::StatementInputTag != data.numberOfPassiveArguments) {
            StmtCallArgs stmtArgs{curAdjointPos, data.numberOfPassiveArguments, curDynamicSizePos, dynamicSizeValues};
            StatementEvaluator::template call<StatementCall::Forward, PrimalValueLinearTape>(
                data.handle, STMT_ARGS_UNPACK(stmtArgs), tape, primalVector, adjointVector, lhsPrimals.data(),
                lhsTangents.data());
          } else {
            curAdjointPos += 1;
          }
        }
      }

      /// \copydoc codi::PrimalValueBaseTape::internalEvaluatePrimal_Step3_EvalStatements
      CODI_INLINE static void internalEvaluatePrimal_Step3_EvalStatements(
          /* data from call */
          PrimalValueLinearTape& tape, Real* primalVector,
          /* data from dynamicSizeData */
          size_t& curDynamicSizePos, size_t const& endDynamicPos, char* const dynamicSizeValues,
          /* data from fixedSizeData */
          size_t& curFixedSizePos, size_t const& endFixedSizePos, char const* const fixedSizeValues,
          /* data from index handler */
          size_t const& startAdjointPos, size_t const& endAdjointPos) {
        CODI_UNUSED(endDynamicPos, endFixedSizePos);

        StackArray<Real> lhsPrimals{};
        StmtFixedDataEntry data;

        size_t curAdjointPos = startAdjointPos;

        while (curAdjointPos < endAdjointPos) {
          curFixedSizePos = data.readForward(fixedSizeValues, curFixedSizePos);

          if (Config::StatementInputTag != data.numberOfPassiveArguments) {
            StmtCallArgs stmtArgs{curAdjointPos, data.numberOfPassiveArguments, curDynamicSizePos, dynamicSizeValues};
            StatementEvaluator::template call<StatementCall::Primal, PrimalValueLinearTape>(
                data.handle, STMT_ARGS_UNPACK(stmtArgs), tape, primalVector, lhsPrimals.data());
          } else {
            curAdjointPos += 1;
          }
        }
      }

      /// \copydoc codi::PrimalValueBaseTape::internalEvaluateReverse_Step3_EvalStatements
      CODI_INLINE static void internalEvaluateReverse_Step3_EvalStatements(
          /* data from call */
          PrimalValueLinearTape& tape, Real* primalVector, ADJOINT_VECTOR_TYPE* adjointVector,
          /* data from dynamicSizeData */
          size_t& curDynamicSizePos, size_t const& endDynamicPos, char* const dynamicSizeValues,
          /* data from fixedSizeData */
          size_t& curFixedSizePos, size_t const& endFixedSizePos, char const* const fixedSizeValues,
          /* data from index handler */
          size_t const& startAdjointPos, size_t const& endAdjointPos) {
        CODI_UNUSED(endDynamicPos, endFixedSizePos);

        StackArray<Gradient> lhsAdjoints{};
        StmtFixedDataEntry data;

        size_t curAdjointPos = startAdjointPos;

        while (curAdjointPos > endAdjointPos) {
          curFixedSizePos = data.readReverse(fixedSizeValues, curFixedSizePos);

          if (Config::StatementInputTag != data.numberOfPassiveArguments) {
            StmtCallArgs stmtArgs{curAdjointPos, data.numberOfPassiveArguments, curDynamicSizePos, dynamicSizeValues};
            StatementEvaluator::template call<StatementCall::Reverse, PrimalValueLinearTape>(
                data.handle, STMT_ARGS_UNPACK(stmtArgs), tape, primalVector, adjointVector, lhsAdjoints.data());
          } else {
            curAdjointPos -= 1;
          }
        }
      }

      /// \copydoc codi::PrimalValueBaseTape::internalResetPrimalValues
      /// Empty implementation; primal values are not overwritten with linear index management.
      CODI_INLINE void internalResetPrimalValues(Position const& pos) {
        CODI_UNUSED(pos);

        // Nothing to do.
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
