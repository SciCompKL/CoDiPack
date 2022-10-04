/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2021 Chair for Scientific Computing (SciComp), TU Kaiserslautern
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
#include <cmath>
#include <type_traits>

#include "../misc/macros.hpp"
#include "../misc/memberStore.hpp"
#include "../config.h"
#include "../expressions/lhsExpressionInterface.hpp"
#include "../expressions/logic/compileTimeTraversalLogic.hpp"
#include "../expressions/logic/helpers/forEachLeafLogic.hpp"
#include "../expressions/logic/helpers/jacobianComputationLogic.hpp"
#include "../expressions/logic/traversalLogic.hpp"
#include "../expressions/referenceActiveType.hpp"
#include "../traits/computationTraits.hpp"
#include "../traits/expressionTraits.hpp"
#include "misc/adjointVectorAccess.hpp"
#include "misc/duplicateJacobianRemover.hpp"
#include "commonTapeImplementation.hpp"
#include "data/chunk.hpp"
#include "data/chunkedData.hpp"
#include "indices/indexManagerInterface.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Type definitions for the Jacobian tapes.
   *
   * @tparam T_Real          See TapeTypesInterface.
   * @tparam T_Gradient      See TapeTypesInterface.
   * @tparam T_IndexManager  Index manager for the tape. Has to implement IndexManagerInterface.
   * @tparam T_Data          See TapeTypesInterface.
   */
  template<typename T_Real, typename T_Gradient, typename T_IndexManager, template<typename, typename> class T_Data>
  struct JacobianTapeTypes : public TapeTypesInterface {
    public:

      using Real = CODI_DD(T_Real, double);                                              ///< See JacobianTapeTypes.
      using Gradient = CODI_DD(T_Gradient, double);                                      ///< See JacobianTapeTypes.
      using IndexManager = CODI_DD(T_IndexManager, CODI_T(IndexManagerInterface<int>));  ///< See JacobianTapeTypes.
      template<typename Chunk, typename Nested>
      using Data = CODI_DD(CODI_T(T_Data<Chunk, Nested>),
                           CODI_T(DataInterface<Nested>));  ///< See JacobianTapeTypes.

      using Identifier = typename IndexManager::Index;  ///< See IndexManagerInterface.

      static bool constexpr IsLinearIndexHandler = IndexManager::IsLinear;  ///< True if the index manager is linear.
      static bool constexpr IsThreadSafeIndexHandler =
            IndexManager::IsThreadSafe;  ///< True if the index manager is thread-safe.
      static bool constexpr IsStaticIndexHandler =
            !IsLinearIndexHandler && !IsThreadSafeIndexHandler; ///< For reuse index management, a static index manager
                                                                ///< is used, unless it is thread-safe.

      /// Statement chunk is either \<argument size\> (linear management) or \<lhs identifier, argument size\>
      /// (reuse management).
      using StatementChunk = typename std::conditional<IsLinearIndexHandler, Chunk1<Config::ArgumentSize>,
                                                       Chunk2<Identifier, Config::ArgumentSize> >::type;
      using StatementData = Data<StatementChunk, IndexManager>;  ///< Statement data vector.

      using JacobianChunk = Chunk2<Real, Identifier>;           ///< Jacobian chunk is \<Jacobian, rhs index\>.
      using JacobianData = Data<JacobianChunk, StatementData>;  ///< Jacobian data vector.

      using NestedData = JacobianData;  ///< See TapeTypesInterface.
  };

  /**
   * @brief Base class for Jacobian tape implementations.
   *
   * This class provides a partial implementation of the FullTapeInterface. Two kinds of functionality have to be added
   * in the implementing classes, namely
   * - methods that manage the vector of adjoint variables and have implications on thread safety,
   * - and methods that depend significantly on the index management scheme and are performance critical.
   *
   * Tape evaluations are performed in three steps with two wrapper steps beforehand. Each methods calls the next
   * method:
   * - evaluate
   * - internalEvaluate*
   * - internalEvaluate*_Step1_ExtFunc
   * - internalEvaluate*_Step2_DataExtraction
   * - internalEvaluate*_Step3_EvalStatements
   * The placeholder stands for Reverse, Forward, or Primal.
   *
   * @tparam T_TapeTypes has to implement JacobianTapeTypes.
   * @tparam T_Impl Type of the final implementation.
   */
  template<typename T_TapeTypes, typename T_Impl>
  struct CommonJacobianTapeImplementation : public CommonTapeImplementation<T_TapeTypes, T_Impl> {
    public:

      /// See JacobianBaseTape.
      using TapeTypes = CODI_DD(
          T_TapeTypes, CODI_T(JacobianTapeTypes<double, double, IndexManagerInterface<int>, DefaultChunkedData>));
      /// See JacobianBaseTape.
      using Impl = CODI_DD(T_Impl, CODI_T(FullTapeInterface<double, double, int, EmptyPosition>));

      using Base = CommonTapeImplementation<TapeTypes, Impl>;  ///< Base class abbreviation.
      friend Base;  ///< Allow the base class to call protected and private methods.

      using Real = typename TapeTypes::Real;                  ///< See TapeTypesInterface.
      using Gradient = typename TapeTypes::Gradient;          ///< See TapeTypesInterface.
      using IndexManager = typename TapeTypes::IndexManager;  ///< See JacobianTapeTypes.
      using Identifier = typename TapeTypes::Identifier;      ///< See TapeTypesInterface.

      using StatementData = typename TapeTypes::StatementData;  ///< See JacobianTapeTypes.
      using JacobianData = typename TapeTypes::JacobianData;    ///< See JacobianTapeTypes.

      using PassiveReal = RealTraits::PassiveReal<Real>;  ///< Basic computation type.

      using NestedPosition = typename JacobianData::Position;  ///< See JacobianTapeTypes.
      using Position = typename Base::Position;                ///< See TapeTypesInterface.

      template<typename Adjoint>
      using VectorAccess =
          AdjointVectorAccess<Real, Identifier, Adjoint>;  ///< Vector access type generated by this tape.

      static bool constexpr AllowJacobianOptimization = true;  ///< See InternalStatementRecordingTapeInterface.
      static bool constexpr HasPrimalValues = false;           ///< See PrimalEvaluationTapeInterface.

      static bool constexpr LinearIndexHandling =
          TapeTypes::IsLinearIndexHandler;                  ///< See IdentifierInformationTapeInterface.
      static bool constexpr RequiresPrimalRestore = false;  ///< See PrimalEvaluationTapeInterface.

    protected:

      using Base::options;

#if CODI_RemoveDuplicateJacobianArguments
      DuplicateJacobianRemover<Real, Identifier> jacobianSorter;  ///< Encapsulates jacobianData to remove duplicated
                                                                  ///< Jacobians.
#endif

      MemberStore<IndexManager, Impl, TapeTypes::IsStaticIndexHandler> indexManager;  ///< Index manager.
      StatementData statementData;  ///< Data stream for statement specific data.
      JacobianData jacobianData;    ///< Data stream for argument specific data.

    private:

      CODI_INLINE Impl const& cast() const {
        return static_cast<Impl const&>(*this);
      }

      CODI_INLINE Impl& cast() {
        return static_cast<Impl&>(*this);
      }

    protected:

      /*******************************************************************************/
      /// @name Interface definition
      /// @{

      /// Perform a forward evaluation of the tape. Arguments are from the recursive eval methods of the DataInterface.
      template<typename... Args>
      static void internalEvaluateForward_Step3_EvalStatements(Args&&... args);

      /// Perform a reverse evaluation of the tape. Arguments are from the recursive eval methods of the DataInterface.
      template<typename... Args>
      static void internalEvaluateReverse_Step3_EvalStatements(Args&&... args);

      /// Add statement specific data to the data streams.
      void pushStmtData(Identifier const& index, Config::ArgumentSize const& numberOfArguments);

      // TODO add the missing adjoint methods

      /// @}

    public:

      /// Constructor
      CommonJacobianTapeImplementation()
          : Base(),
#if CODI_RemoveDuplicateJacobianArguments
            jacobianSorter(),
#endif
            indexManager(0),  // Reserve the zero index.
            statementData(Config::ChunkSize),
            jacobianData(Config::ChunkSize)
      {
        statementData.setNested(&indexManager.get());
        jacobianData.setNested(&statementData);

        Base::init(&jacobianData);

        Base::options.insert(TapeParameters::JacobianSize);
        Base::options.insert(TapeParameters::LargestIdentifier);
        Base::options.insert(TapeParameters::StatementSize);
      }

      /*******************************************************************************/
      /// @name Functions from InternalStatementRecordingTapeInterface
      /// @{

      /// \copydoc codi::InternalStatementRecordingTapeInterface::initIdentifier()
      template<typename Real>
      CODI_INLINE void initIdentifier(Real& value, Identifier& identifier) {
        CODI_UNUSED(value);

        identifier = IndexManager::InactiveIndex;
      }

      /// \copydoc codi::InternalStatementRecordingTapeInterface::destroyIdentifier()
      template<typename Real>
      CODI_INLINE void destroyIdentifier(Real& value, Identifier& identifier) {
        CODI_UNUSED(value);

        indexManager.get().freeIndex(identifier);
      }

      /// @}

    protected:

      /// Pushes Jacobians and indices to the tape.
      struct PushJacobianLogic : public JacobianComputationLogic<PushJacobianLogic> {
        public:
          /// General implementation. Checks for invalid and passive values/Jacobians.
          template<typename Node, typename Jacobian, typename DataVector>
          CODI_INLINE void handleJacobianOnActive(Node const& node, Jacobian jacobianExpr, DataVector& dataVector) {
            Real jacobian = ComputationTraits::adjointConversion<Real>(jacobianExpr);

            if (CODI_ENABLE_CHECK(Config::CheckZeroIndex, 0 != node.getIdentifier())) {
              if (CODI_ENABLE_CHECK(Config::IgnoreInvalidJacobians, RealTraits::isTotalFinite(jacobian))) {
                if (CODI_ENABLE_CHECK(Config::CheckJacobianIsZero, !RealTraits::isTotalZero(jacobian))) {
                  dataVector.pushData(jacobian, node.getIdentifier());
                }
              }
            }
          }

          /// Specialization for ReferenceActiveType nodes. Delays Jacobian push.
          template<typename Type, typename Jacobian, typename DataVector>
          CODI_INLINE void handleJacobianOnActive(ReferenceActiveType<Type> const& node, Jacobian jacobianExpr,
                                                  DataVector& dataVector) {
            CODI_UNUSED(dataVector);

            Real jacobian = ComputationTraits::adjointConversion<Real>(jacobianExpr);

            if (CODI_ENABLE_CHECK(Config::IgnoreInvalidJacobians, RealTraits::isTotalFinite(jacobian))) {
              // Do a delayed push for these leaf nodes, accumulate the jacobian in the local member.
              node.jacobian += jacobian;
            }
          }
      };

      /// Pushes all delayed Jacobians.
      struct PushDelayedJacobianLogic : public ForEachLeafLogic<PushDelayedJacobianLogic> {
        public:

          /// Specialization for ReferenceActiveType nodes. Pushes the delayed Jacobian.
          template<typename Type, typename DataVector>
          CODI_INLINE void handleActive(ReferenceActiveType<Type> const& node, DataVector& dataVector) {
            if (CODI_ENABLE_CHECK(Config::CheckZeroIndex, 0 != node.getIdentifier())) {
              if (CODI_ENABLE_CHECK(Config::CheckJacobianIsZero, !RealTraits::isTotalZero(node.jacobian))) {
                dataVector.pushData(node.jacobian, node.getIdentifier());

                // Reset the jacobian here such that it is not pushed multiple times and ready for the next store.
                node.jacobian = Real();
              }
            }
          }

          using ForEachLeafLogic<PushDelayedJacobianLogic>::handleActive;
      };

      /// Push Jacobians and delayed Jacobians to the tape.
      template<typename Rhs>
      CODI_INLINE void pushJacobians(ExpressionInterface<Real, Rhs> const& rhs) {
        PushJacobianLogic pushJacobianLogic;
        PushDelayedJacobianLogic pushDelayedJacobianLogic;

#if CODI_RemoveDuplicateJacobianArguments
        auto& insertVector = jacobianSorter;
#else
        auto& insertVector = jacobianData;
#endif

        pushJacobianLogic.eval(rhs.cast(), Real(1.0), insertVector);
        pushDelayedJacobianLogic.eval(rhs.cast(), insertVector);

#if CODI_RemoveDuplicateJacobianArguments
        jacobianSorter.storeData(jacobianData);
#endif
      }

    public:

      /// @{

      /// \copydoc codi::InternalStatementRecordingTapeInterface::store()
      template<typename Lhs, typename Rhs>
      CODI_INLINE void store(LhsExpressionInterface<Real, Gradient, Impl, Lhs>& lhs,
                             ExpressionInterface<Real, Rhs> const& rhs) {
        if (CODI_ENABLE_CHECK(Config::CheckTapeActivity, cast().isActive())) {
          size_t constexpr MaxArgs = ExpressionTraits::NumberOfActiveTypeArguments<Rhs>::value;

          codiAssert(MaxArgs < Config::MaxArgumentSize);

          statementData.reserveItems(1);
          typename JacobianData::InternalPosHandle jacobianStart = jacobianData.reserveItems(MaxArgs);

          pushJacobians(rhs);

          size_t numberOfArguments = jacobianData.getPushedDataCount(jacobianStart);
          if (CODI_ENABLE_CHECK(Config::CheckEmptyStatements, 0 != numberOfArguments)) {
            indexManager.get().assignIndex(lhs.cast().getIdentifier());
            cast().pushStmtData(lhs.cast().getIdentifier(), (Config::ArgumentSize)numberOfArguments);
          } else {
            indexManager.get().freeIndex(lhs.cast().getIdentifier());
          }
        } else {
          indexManager.get().freeIndex(lhs.cast().getIdentifier());
        }

        lhs.cast().value() = rhs.cast().getValue();
      }

      /// \copydoc codi::InternalStatementRecordingTapeInterface::store() <br>
      /// Optimization for copy statements.
      template<typename Lhs, typename Rhs>
      CODI_INLINE void store(LhsExpressionInterface<Real, Gradient, Impl, Lhs>& lhs,
                             LhsExpressionInterface<Real, Gradient, Impl, Rhs> const& rhs) {
        if (CODI_ENABLE_CHECK(Config::CheckTapeActivity, cast().isActive())) {
          if (IndexManager::CopyNeedsStatement || !Config::CopyOptimization) {
            store<Lhs, Rhs>(lhs, static_cast<ExpressionInterface<Real, Rhs> const&>(rhs));
          } else {
            indexManager.get().copyIndex(lhs.cast().getIdentifier(), rhs.cast().getIdentifier());
          }
        } else {
          indexManager.get().freeIndex(lhs.cast().getIdentifier());
        }

        lhs.cast().value() = rhs.cast().getValue();
      }

      /// \copydoc codi::InternalStatementRecordingTapeInterface::store() <br>
      /// Specialization for passive assignments.
      template<typename Lhs>
      CODI_INLINE void store(LhsExpressionInterface<Real, Gradient, Impl, Lhs>& lhs, Real const& rhs) {
        indexManager.get().freeIndex(lhs.cast().getIdentifier());

        lhs.cast().value() = rhs;
      }

      /// @}
      /*******************************************************************************/
      /// @name Functions from ReverseTapeInterface
      /// @{

      /// \copydoc codi::ReverseTapeInterface::registerInput()
      template<typename Lhs>
      CODI_INLINE void registerInput(LhsExpressionInterface<Real, Gradient, Impl, Lhs>& value) {
        indexManager.get().assignUnusedIndex(value.cast().getIdentifier());

        if (TapeTypes::IsLinearIndexHandler) {
          statementData.reserveItems(1);
          cast().pushStmtData(value.cast().getIdentifier(), Config::StatementInputTag);
        }
      }

      /// @}

    protected:

      /// Adds data from the index, statement, and Jacobian streams to the tape values.
      CODI_INLINE void internalAddTapeValues(TapeValues& values) const {
        values.addSection("Index manager");
        indexManager.get().addToTapeValues(values);
        values.addSection("Statement entries");
        statementData.addToTapeValues(values);
        values.addSection("Jacobian entries");
        jacobianData.addToTapeValues(values);
      }

      /******************************************************************************
       * Protected helper function for CustomAdjointVectorEvaluationTapeInterface
       */

    protected:

      /// Performs the AD \ref sec_reverseAD "reverse" equation for a statement.
      template<typename Adjoint>
      CODI_INLINE static void incrementAdjoints(Adjoint* adjointVector, Adjoint const& lhsAdjoint,
                                                Config::ArgumentSize const& numberOfArguments, size_t& curJacobianPos,
                                                Real const* const rhsJacobians,
                                                Identifier const* const rhsIdentifiers) {
        size_t endJacobianPos = curJacobianPos - numberOfArguments;

        if (CODI_ENABLE_CHECK(Config::SkipZeroAdjointEvaluation, !RealTraits::isTotalZero(lhsAdjoint))) {
          while (endJacobianPos < curJacobianPos) {
            curJacobianPos -= 1;
            adjointVector[rhsIdentifiers[curJacobianPos]] += rhsJacobians[curJacobianPos] * lhsAdjoint;
          }
        } else {
          curJacobianPos = endJacobianPos;
        }
      }

      /// Wrapper helper for improved compiler optimizations.
      CODI_WRAP_FUNCTION_TEMPLATE(Wrap_internalEvaluateReverse_Step3_EvalStatements,
                                  Impl::template internalEvaluateReverse_Step3_EvalStatements);

      /// Start for reverse evaluation between external function.
      template<typename Adjoint>
      CODI_NO_INLINE static void internalEvaluateReverse_Step2_DataExtraction(NestedPosition const& start,
                                                                              NestedPosition const& end, Adjoint* data,
                                                                              JacobianData& jacobianData) {
        Wrap_internalEvaluateReverse_Step3_EvalStatements<Adjoint> evalFunc;
        jacobianData.evaluateReverse(start, end, evalFunc, data);
      }

      /// Performs the AD \ref sec_forwardAD "forward" equation for a statement.
      template<typename Adjoint>
      CODI_INLINE static void incrementTangents(Adjoint const* const adjointVector, Adjoint& lhsAdjoint,
                                                Config::ArgumentSize const& numberOfArguments, size_t& curJacobianPos,
                                                Real const* const rhsJacobians,
                                                Identifier const* const rhsIdentifiers) {
        size_t endJacobianPos = curJacobianPos + numberOfArguments;

        while (curJacobianPos < endJacobianPos) {
          lhsAdjoint += rhsJacobians[curJacobianPos] * adjointVector[rhsIdentifiers[curJacobianPos]];
          curJacobianPos += 1;
        }
      }

      /// Wrapper helper for improved compiler optimizations.
      CODI_WRAP_FUNCTION_TEMPLATE(Wrap_internalEvaluateForward_Step3_EvalStatements,
                                  Impl::template internalEvaluateForward_Step3_EvalStatements);

      /// Start for forward evaluation between external function.
      template<typename Adjoint>
      CODI_NO_INLINE static void internalEvaluateForward_Step2_DataExtraction(NestedPosition const& start,
                                                                              NestedPosition const& end, Adjoint* data,
                                                                              JacobianData& jacobianData) {
        Wrap_internalEvaluateForward_Step3_EvalStatements<Adjoint> evalFunc;
        jacobianData.evaluateForward(start, end, evalFunc, data);
      }

    public:

      /// @name Functions from CustomAdjointVectorEvaluationTapeInterface
      /// @{

      using Base::evaluate;

      /// \copydoc codi::CustomAdjointVectorEvaluationTapeInterface::evaluate()
      template<typename Adjoint>
      CODI_NO_INLINE void evaluate(Position const& start, Position const& end, Adjoint* data) {
        VectorAccess<Adjoint> adjointWrapper(data);

        Base::internalEvaluateReverse_Step1_ExtFunc(
            start, end,
            CommonJacobianTapeImplementation::template internalEvaluateReverse_Step2_DataExtraction<Adjoint>,
            &adjointWrapper, data, jacobianData);
      }

      /// \copydoc codi::CustomAdjointVectorEvaluationTapeInterface::evaluate()
      template<typename Adjoint>
      CODI_NO_INLINE void evaluateForward(Position const& start, Position const& end, Adjoint* data) {
        VectorAccess<Adjoint> adjointWrapper(data);

        Base::internalEvaluateForward_Step1_ExtFunc(
            start, end,
            CommonJacobianTapeImplementation::template internalEvaluateForward_Step2_DataExtraction<Adjoint>,
            &adjointWrapper, data, jacobianData);
      }

      /// @}
      /*******************************************************************************/
      /// @name Functions from DataManagementTapeInterface
      /// @{

      /// \copydoc codi::DataManagementTapeInterface::swap()
      CODI_INLINE void swap(Impl& other) {
        // Index manager does not need to be swapped, it is either static or swapped with the vector data.
        // Vectors are swapped recursively in the base class.

        Base::swap(other);
      }

      /// \copydoc codi::DataManagementTapeInterface::getParameter()
      size_t getParameter(TapeParameters parameter) const {
        switch (parameter) {
          case TapeParameters::JacobianSize:
            return jacobianData.getDataSize();
            break;
          case TapeParameters::LargestIdentifier:
            return indexManager.get().getLargestCreatedIndex();
            break;
          case TapeParameters::StatementSize:
            return statementData.getDataSize();
          default:
            return Base::getParameter(parameter);
            break;
        }
      }

      /// \copydoc codi::DataManagementTapeInterface::setParameter()
      void setParameter(TapeParameters parameter, size_t value) {
        switch (parameter) {
          case TapeParameters::JacobianSize:
            jacobianData.resize(value);
            break;
          case TapeParameters::LargestIdentifier:
            CODI_EXCEPTION("Tried to set a get only parameter.");
            break;
          case TapeParameters::StatementSize:
            statementData.resize(value);
            break;
          default:
            Base::setParameter(parameter, value);
            break;
        }
      }

      /// @}
      /*******************************************************************************/
      /// @name Functions from ExternalFunctionTapeInterface
      /// @{

      /// \copydoc codi::ExternalFunctionTapeInterface::registerExternalFunctionOutput()
      template<typename Lhs>
      Real registerExternalFunctionOutput(LhsExpressionInterface<Real, Gradient, Impl, Lhs>& value) {
        cast().registerInput(value);

        return Real();
      }

      /// @}
      /*******************************************************************************/
      /// @name Functions from ManualStatementPushTapeInterface
      /// @{

      /// \copydoc codi::ManualStatementPushTapeInterface::pushJacobiManual()
      void pushJacobiManual(Real const& jacobian, Real const& value, Identifier const& index) {
        CODI_UNUSED(value);

        jacobianData.pushData(jacobian, index);
      }

      /// \copydoc codi::ManualStatementPushTapeInterface::storeManual()
      void storeManual(Real const& lhsValue, Identifier& lhsIndex, Config::ArgumentSize const& size) {
        CODI_UNUSED(lhsValue);

        codiAssert(size < Config::MaxArgumentSize);

        statementData.reserveItems(1);
        jacobianData.reserveItems(size);

        indexManager.get().assignIndex(lhsIndex);
        cast().pushStmtData(lhsIndex, (Config::ArgumentSize)size);
      }

      /// @}
      /*******************************************************************************/
      /// @name Functions from PreaccumulationEvaluationTapeInterface
      /// @{

      /// \copydoc codi::PreaccumulationEvaluationTapeInterface::evaluateKeepState()
      void evaluateKeepState(Position const& start, Position const& end) {
        cast().evaluate(start, end);
      }

      /// \copydoc codi::PreaccumulationEvaluationTapeInterface::evaluateForwardKeepState()
      void evaluateForwardKeepState(Position const& start, Position const& end) {
        cast().evaluateForward(start, end);
      }

      /// @}
      /*******************************************************************************/
      /// @name Functions from PrimalEvaluationTapeInterface
      /// @{

      using Base::evaluatePrimal;

      /// Not implemented, raises an exception.
      void evaluatePrimal(Position const& start, Position const& end) {
        CODI_UNUSED(start, end);

        CODI_EXCEPTION("Accessing primal evaluation of an Jacobian tape.");
      }

      /// Not implemented, raises an exception.
      Real& primal(Identifier const& identifier) {
        CODI_UNUSED(identifier);

        CODI_EXCEPTION("Accessing primal vector of an Jacobian tape.");

        static Real temp;
        return temp;
      }

      /// Not implemented, raises an exception.
      Real primal(Identifier const& identifier) const {
        CODI_UNUSED(identifier);

        CODI_EXCEPTION("Accessing primal vector of an Jacobian tape.");

        return Real();
      }

      /// @}

  };
}
