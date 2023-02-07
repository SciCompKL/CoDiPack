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

      using StmtFixedDataEntry = typename Base::StmtFixedDataEntry;  ///< See PrimalValueBaseTape.
      using StmtCallArgs = typename Base::StmtCallArgs;              ///< See PrimalValueBaseTape.

      template<typename T>
      using StackArray = typename Base::template StackArray<T>;  ///< See PrimalValueBaseTape.

      /// Constructor
      PrimalValueReuseTape() : Base() {}

      using Base::clearAdjoints;

      /// \copydoc codi::PositionalEvaluationTapeInterface::clearAdjoints
      void clearAdjoints(Position const& start, Position const& end) {
        auto clearFunc = [](
                             /* data from call */
                             PrimalValueReuseTape &tape, ADJOINT_VECTOR_TYPE* adjointVector, size_t adjointVectorSize,
                             /* data from dynamicSizeData */
                             size_t& curDynamicSizePos, size_t const& endDynamicPos, char* const dynamicSizeValues,
                             /* data from fixedSizeData */
                             size_t& curFixedSizePos, size_t const& endFixedSizePos,
                             char const* const fixedSizeValues) {
          CODI_UNUSED(endDynamicPos);

          StmtFixedDataEntry data;

          while (curFixedSizePos > endFixedSizePos) {
            curFixedSizePos = data.readReverse(fixedSizeValues, curFixedSizePos);

            StmtCallArgs stmtArgs{data.numberOfPassiveArguments, curDynamicSizePos, dynamicSizeValues};
            StatementEvaluator::template call<StatementCall::ClearAdjoint, PrimalValueReuseTape>(
                data.handle, STMT_ARGS_UNPACK(stmtArgs), tape, adjointVector, adjointVectorSize);
          }
        };

        using DynamicPosition = typename Base::DynamicSizeData::Position;
        DynamicPosition startStmt = this->externalFunctionData.template extractPosition<DynamicPosition>(start);
        DynamicPosition endStmt = this->externalFunctionData.template extractPosition<DynamicPosition>(end);

        typename Base::template VectorAccess<Gradient> vectorAccess(Base::adjoints.data(), Base::primals.data());

        ADJOINT_VECTOR_TYPE* dataVector = Base::selectAdjointVector(&vectorAccess, Base::adjoints.data());

        this->dynamicSizeData.evaluateReverse(startStmt, endStmt, clearFunc, *this, dataVector, Base::adjoints.size());
      }

    protected:

      /// \copydoc codi::PrimalValueBaseTape::internalEvaluateForward_Step3_EvalStatements
      CODI_INLINE static void internalEvaluateForward_Step3_EvalStatements(
          /* data from call */
          PrimalValueReuseTape& tape, Real* primalVector, ADJOINT_VECTOR_TYPE* adjointVector,
          /* data from dynamicSizeData */
          size_t& curDynamicSizePos, size_t const& endDynamicPos, char* const dynamicSizeValues,
          /* data from fixedSizeData */
          size_t& curFixedSizePos, size_t const& endFixedSizePos, char const* const fixedSizeValues) {
        CODI_UNUSED(endDynamicPos);

        StackArray<Real> lhsPrimals{};
        StackArray<Gradient> lhsTangents{};
        StmtFixedDataEntry data;

        while (curFixedSizePos < endFixedSizePos) {
          curFixedSizePos = data.readForward(fixedSizeValues, curFixedSizePos);

          StmtCallArgs stmtArgs{data.numberOfPassiveArguments, curDynamicSizePos, dynamicSizeValues};
          StatementEvaluator::template call<StatementCall::Forward, PrimalValueReuseTape>(
              data.handle, STMT_ARGS_UNPACK(stmtArgs), tape, primalVector, adjointVector, lhsPrimals.data(),
              lhsTangents.data());
        }
      }

      /// \copydoc codi::PrimalValueBaseTape::internalEvaluatePrimal_Step3_EvalStatements
      CODI_INLINE static void internalEvaluatePrimal_Step3_EvalStatements(
          /* data from call */
          PrimalValueReuseTape& tape, Real* primalVector,
          /* data from dynamicSizeData */
          size_t& curDynamicSizePos, size_t const& endDynamicPos, char* const dynamicSizeValues,
          /* data from fixedSizeData */
          size_t& curFixedSizePos, size_t const& endFixedSizePos, char const* const fixedSizeValues) {
        CODI_UNUSED(endDynamicPos);

        StackArray<Real> lhsPrimals{};
        StmtFixedDataEntry data;

        while (curFixedSizePos < endFixedSizePos) {
          curFixedSizePos = data.readForward(fixedSizeValues, curFixedSizePos);

          StmtCallArgs stmtArgs{data.numberOfPassiveArguments, curDynamicSizePos, dynamicSizeValues};
          StatementEvaluator::template call<StatementCall::Primal, PrimalValueReuseTape>(
              data.handle, STMT_ARGS_UNPACK(stmtArgs), tape, primalVector, lhsPrimals.data());
        }
      }

      /// \copydoc codi::PrimalValueBaseTape::internalEvaluateReverse_Step3_EvalStatements
      CODI_NO_INLINE static void internalEvaluateReverse_Step3_EvalStatements(
          /* data from call */
          PrimalValueReuseTape& tape, Real* primalVector, ADJOINT_VECTOR_TYPE* adjointVector,
          /* data from dynamicSizeData */
          size_t& curDynamicSizePos, size_t const& endDynamicPos, char* const dynamicSizeValues,
          /* data from fixedSizeData */
          size_t& curFixedSizePos, size_t const& endFixedSizePos, char const* const fixedSizeValues) {
        CODI_UNUSED(endDynamicPos);

        StackArray<Gradient> lhsAdjoints{};
        StmtFixedDataEntry data;

        while (curFixedSizePos > endFixedSizePos) {
          curFixedSizePos = data.readReverse(fixedSizeValues, curFixedSizePos);

          StmtCallArgs stmtArgs{data.numberOfPassiveArguments, curDynamicSizePos, dynamicSizeValues};
          StatementEvaluator::template call<StatementCall::Reverse, PrimalValueReuseTape>(
              data.handle, STMT_ARGS_UNPACK(stmtArgs), tape, primalVector, adjointVector, lhsAdjoints.data());
        }
      }

      /// \copydoc codi::PrimalValueBaseTape::internalResetPrimalValues
      CODI_INLINE void internalResetPrimalValues(Position const& pos) {
        auto clearFunc = [](
                             /* data from call */
                             PrimalValueReuseTape &tape, Real* primalVector,
                             /* data from dynamicSizeData */
                             size_t& curDynamicSizePos, size_t const& endDynamicPos, char* const dynamicSizeValues,
                             /* data from fixedSizeData */
                             size_t& curFixedSizePos, size_t const& endFixedSizePos,
                             char const* const fixedSizeValues) {
          CODI_UNUSED(endDynamicPos);

          StmtFixedDataEntry data;

          while (curFixedSizePos > endFixedSizePos) {
            curFixedSizePos = data.readReverse(fixedSizeValues, curFixedSizePos);

            StmtCallArgs stmtArgs{data.numberOfPassiveArguments, curDynamicSizePos, dynamicSizeValues};
            StatementEvaluator::template call<StatementCall::ResetPrimal, PrimalValueReuseTape>(
                data.handle, STMT_ARGS_UNPACK(stmtArgs), tape, primalVector);
          }
        };

        using DynamicPosition = typename Base::DynamicSizeData::Position;
        DynamicPosition startStmt =
            this->externalFunctionData.template extractPosition<DynamicPosition>(this->getPosition());
        DynamicPosition endStmt = this->externalFunctionData.template extractPosition<DynamicPosition>(pos);

        this->dynamicSizeData.evaluateReverse(startStmt, endStmt, clearFunc, *this, this->primals.data());
      }

    public:
      /// \copydoc codi::PrimalEvaluationTapeInterface::revertPrimals
      void revertPrimals(Position const& pos) {
        internalResetPrimalValues(pos);
      }
  };
}
