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
#include <cmath>
#include <type_traits>

#include "../config.h"
#include "../expressions/lhsExpressionInterface.hpp"
#include "../expressions/logic/compileTimeTraversalLogic.hpp"
#include "../expressions/logic/helpers/forEachLeafLogic.hpp"
#include "../expressions/logic/helpers/jacobianComputationLogic.hpp"
#include "../expressions/logic/traversalLogic.hpp"
#include "../expressions/referenceActiveType.hpp"
#include "../misc/macros.hpp"
#include "../misc/mathUtility.hpp"
#include "../misc/memberStore.hpp"
#include "../traits/computationTraits.hpp"
#include "../traits/expressionTraits.hpp"
#include "commonTapeImplementation.hpp"
#include "data/chunk.hpp"
#include "data/chunkedData.hpp"
#include "indices/indexManagerInterface.hpp"
#include "misc/adjointVectorAccess.hpp"
#include "misc/duplicateJacobianRemover.hpp"
#include "misc/localAdjoints.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Type definitions for the Jacobian tapes.
   *
   * @tparam T_Real          See TapeTypesInterface.
   * @tparam T_Gradient      See TapeTypesInterface.
   * @tparam T_IndexManager  Index manager for the tape. Has to implement IndexManagerInterface.
   * @tparam T_Data          See TapeTypesInterface.
   * @tparam T_Adjoints      Internal implementation of the adjoint variables.
   */
  template<typename T_Real, typename T_Gradient, typename T_IndexManager, template<typename, typename> class T_Data,
           template<typename, typename, typename> class T_Adjoints = LocalAdjoints>
  struct JacobianTapeTypes : public TapeTypesInterface {
    public:

      using Real = CODI_DD(T_Real, double);                                              ///< See JacobianTapeTypes.
      using Gradient = CODI_DD(T_Gradient, double);                                      ///< See JacobianTapeTypes.
      using IndexManager = CODI_DD(T_IndexManager, CODI_T(IndexManagerInterface<int>));  ///< See JacobianTapeTypes.
      template<typename Chunk, typename Nested>
      using Data = CODI_DD(CODI_T(T_Data<Chunk, Nested>),
                           CODI_T(DataInterface<Nested>));  ///< See JacobianTapeTypes.

      using Identifier = typename IndexManager::Index;  ///< See IndexManagerInterface.

      /// See JacobianTapeTypes.
      template<typename Impl>
      using Adjoints = CODI_DD(CODI_T(T_Adjoints<Gradient, Identifier, Impl>),
                               CODI_T(InternalAdjointsInterface<Gradient, Identifier, CODI_ANY>));

      static bool constexpr IsLinearIndexHandler = IndexManager::IsLinear;  ///< True if the index manager is linear.
      static bool constexpr IsStaticIndexHandler =
          IndexManager::NeedsStaticStorage;  ///< True if the index manager must be stored statically in the tape.

      /// Statement chunk is either \<argument size\> (linear management) or \<lhs identifier, argument size\>
      /// (reuse management).
      using StatementChunk = typename std::conditional<IsLinearIndexHandler, Chunk1<Config::ArgumentSize>,
                                                       Chunk2<Identifier, Config::ArgumentSize> >::type;
      using StatementData = Data<StatementChunk, IndexManager>;  ///< Statement data vector.

      using JacobianChunk = Chunk2<Real, Identifier>;           ///< Jacobian chunks is \<Jacobian, rhs index\>.
      using JacobianData = Data<JacobianChunk, StatementData>;  ///< Jacobian data vector.

      using NestedData = JacobianData;  ///< See TapeTypesInterface.
  };

  /**
   * @brief Base class for all standard Jacobian tape implementations.
   *
   * This class provides nearly a full implementation of the FullTapeInterface. There are just a few internal methods
   * left which have to be implemented by the final classes. These methods depend significantly on the index management
   * scheme and are performance critical.
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
  struct JacobianBaseTape : public CommonTapeImplementation<T_TapeTypes, T_Impl> {
    public:

      /// See JacobianBaseTape.
      using TapeTypes = CODI_DD(
          T_TapeTypes, CODI_T(JacobianTapeTypes<double, double, IndexManagerInterface<int>, DefaultChunkedData>));
      using Impl = CODI_DD(T_Impl, CODI_DEFAULT_TAPE);  ///< See JacobianBaseTape.

      using Base = CommonTapeImplementation<T_TapeTypes, T_Impl>;  ///< Base class abbreviation.
      friend Base;  ///< Allow the base class to call protected and private methods.

      using Real = typename TapeTypes::Real;                  ///< See TapeTypesInterface.
      using Gradient = typename TapeTypes::Gradient;          ///< See TapeTypesInterface.
      using IndexManager = typename TapeTypes::IndexManager;  ///< See JacobianTapeTypes.
      using Identifier = typename TapeTypes::Identifier;      ///< See TapeTypesInterface.

      using StatementData = typename TapeTypes::StatementData;  ///< See JacobianTapeTypes.
      using JacobianData = typename TapeTypes::JacobianData;    ///< See JacobianTapeTypes.

      using Adjoints = typename TapeTypes::template Adjoints<Impl>;  ///< See JacobianTapeTypes.

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

#if CODI_RemoveDuplicateJacobianArguments
      DuplicateJacobianRemover<Real, Identifier> jacobianSorter;  ///< Encapsulates jacobianData to remove duplicated
                                                                  ///< Jacobians.
#endif

      MemberStore<IndexManager, Impl, TapeTypes::IsStaticIndexHandler> indexManager;  ///< Index manager.
      StatementData statementData;  ///< Data stream for statement specific data.
      JacobianData jacobianData;    ///< Data stream for argument specific data.

      Adjoints adjoints;  ///< Evaluation vector for AD.

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

      /// @}

    public:

      /// Constructor
      JacobianBaseTape()
          : Base(),
#if CODI_RemoveDuplicateJacobianArguments
            jacobianSorter(),
#endif
            indexManager(0),  // Reserve the zero index.
            statementData(Config::ChunkSize),
            jacobianData(std::max(Config::ChunkSize, Config::MaxArgumentSize)),  // Chunk must be large enough to store
                                                                                 // data for all arguments of one
                                                                                 // statement.
            adjoints(1)  // Ensure that adjoint[0] exists, see its use in gradient() const.
      {
        statementData.setNested(&indexManager.get());
        jacobianData.setNested(&statementData);

        Base::init(&jacobianData);

        Base::options.insert(TapeParameters::AdjointSize);
        Base::options.insert(TapeParameters::JacobianSize);
        Base::options.insert(TapeParameters::LargestIdentifier);
        Base::options.insert(TapeParameters::StatementSize);
      }

      /*******************************************************************************/
      /// @name Functions from GradientAccessTapeInterface
      /// @{

      /// \copydoc codi::GradientAccessTapeInterface::gradient(Identifier const&, AdjointsManagement)
      CODI_INLINE Gradient& gradient(Identifier const& identifier,
                                     AdjointsManagement adjointsManagement = AdjointsManagement::Automatic) {
        if (AdjointsManagement::Automatic == adjointsManagement) {
          checkAdjointSize(identifier);
        }

        codiAssert(identifier < (Identifier)adjoints.size());

        return adjoints[identifier];
      }

      /// \copydoc codi::GradientAccessTapeInterface::gradient(Identifier const&, AdjointsManagement) const
      CODI_INLINE Gradient const& gradient(
          Identifier const& identifier, AdjointsManagement adjointsManagement = AdjointsManagement::Automatic) const {
        codiAssert(identifier < (Identifier)adjoints.size());

        if (AdjointsManagement::Automatic == adjointsManagement && identifier >= (Identifier)adjoints.size()) {
          return adjoints[0];
        } else {
          return adjoints[identifier];
        }
      }

      /// @}
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

        indexManager.get().template freeIndex<Impl>(identifier);
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
            indexManager.get().template assignIndex<Impl>(lhs.cast().getIdentifier());
            cast().pushStmtData(lhs.cast().getIdentifier(), (Config::ArgumentSize)numberOfArguments);

            if (Config::StatementEvents) {
              Real* jacobians;
              Identifier* rhsIdentifiers;
              jacobianData.getDataPointers(jacobianStart, jacobians, rhsIdentifiers);

              EventSystem<Impl>::notifyStatementStoreOnTapeListeners(cast(), lhs.cast().getIdentifier(),
                                                                     rhs.cast().getValue(), numberOfArguments,
                                                                     rhsIdentifiers, jacobians);
            }
          } else {
            indexManager.get().template freeIndex<Impl>(lhs.cast().getIdentifier());
          }
        } else {
          indexManager.get().template freeIndex<Impl>(lhs.cast().getIdentifier());
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
            return;
          } else {
            indexManager.get().template copyIndex<Impl>(lhs.cast().getIdentifier(), rhs.cast().getIdentifier());
          }
        } else {
          indexManager.get().template freeIndex<Impl>(lhs.cast().getIdentifier());
        }

        lhs.cast().value() = rhs.cast().getValue();
      }

      /// \copydoc codi::InternalStatementRecordingTapeInterface::store() <br>
      /// Specialization for passive assignments.
      template<typename Lhs>
      CODI_INLINE void store(LhsExpressionInterface<Real, Gradient, Impl, Lhs>& lhs, Real const& rhs) {
        indexManager.get().template freeIndex<Impl>(lhs.cast().getIdentifier());

        lhs.cast().value() = rhs;
      }

      /// @}
      /*******************************************************************************
       * Protected helper function for ReverseTapeInterface
       */

    protected:

      /// Add a new input to the tape.
      template<typename Lhs>
      CODI_INLINE void internalRegisterInput(LhsExpressionInterface<Real, Gradient, Impl, Lhs>& value,
                                             bool unusedIndex) {
        if (TapeTypes::IsLinearIndexHandler) {
          statementData.reserveItems(1);
        }

        if (unusedIndex) {
          indexManager.get().template assignUnusedIndex<Impl>(value.cast().getIdentifier());
        } else {
          indexManager.get().template assignIndex<Impl>(value.cast().getIdentifier());
        }

        if (TapeTypes::IsLinearIndexHandler) {
          cast().pushStmtData(value.cast().getIdentifier(), Config::StatementInputTag);
        }
      }

    public:

      /*******************************************************************************/
      /// @name Functions from ReverseTapeInterface
      /// @{

      /// \copydoc codi::ReverseTapeInterface::registerInput()
      template<typename Lhs>
      CODI_INLINE void registerInput(LhsExpressionInterface<Real, Gradient, Impl, Lhs>& value) {
        cast().internalRegisterInput(value, true);
        EventSystem<Impl>::notifyTapeRegisterInputListeners(cast(), value.cast().value(), value.cast().getIdentifier());
      }

      /// \copydoc codi::ReverseTapeInterface::clearAdjoints()
      CODI_INLINE void clearAdjoints(AdjointsManagement adjointsManagement = AdjointsManagement::Automatic) {
        if (AdjointsManagement::Automatic == adjointsManagement) {
          adjoints.beginUse();
        }

        adjoints.zeroAll();

        if (AdjointsManagement::Automatic == adjointsManagement) {
          adjoints.endUse();
        }
      }

      /// @}

    protected:

      /// Adds data from all streams, the size of the adjoint vector and index manager data.
      CODI_INLINE TapeValues internalGetTapeValues() const {
        std::string name;
        if (TapeTypes::IsLinearIndexHandler) {
          name = "CoDi Tape Statistics ( JacobianLinearTape )";
        } else {
          name = "CoDi Tape Statistics ( JacobianReuseTape )";
        }
        TapeValues values = TapeValues(name);

        size_t nAdjoints = indexManager.get().getLargestCreatedIndex();
        double memoryAdjoints = static_cast<double>(nAdjoints) * static_cast<double>(sizeof(Gradient));

        values.addSection("Adjoint vector");
        values.addUnsignedLongEntry("Number of adjoints", nAdjoints);
        values.addDoubleEntry("Memory allocated", memoryAdjoints, true, true);

        values.addSection("Index manager");
        indexManager.get().addToTapeValues(values);

        values.addSection("Statement entries");
        statementData.addToTapeValues(values);
        values.addSection("Jacobian entries");
        jacobianData.addToTapeValues(values);

        return values;
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

        if (CODI_ENABLE_CHECK(Config::SkipZeroAdjointEvaluation, !RealTraits::isTotalZero(lhsAdjoint))) CODI_Likely {
          while (endJacobianPos < curJacobianPos) CODI_Likely {
            curJacobianPos -= 1;
            adjointVector[rhsIdentifiers[curJacobianPos]] += rhsJacobians[curJacobianPos] * lhsAdjoint;
          }
        } else CODI_Unlikely {
          curJacobianPos = endJacobianPos;
        }
      }

      /// Wrapper helper for improved compiler optimizations.
      CODI_WRAP_FUNCTION_TEMPLATE(Wrap_internalEvaluateReverse_Step3_EvalStatements,
                                  Impl::template internalEvaluateReverse_Step3_EvalStatements);

      /// Start for reverse evaluation between external function.
      template<typename Adjoint>
      CODI_NO_INLINE static void internalEvaluateReverse_Step2_DataExtraction(NestedPosition const& start,
                                                                              NestedPosition const& end, Impl& tape,
                                                                              Adjoint* data,
                                                                              JacobianData& jacobianData) {
        Wrap_internalEvaluateReverse_Step3_EvalStatements<Adjoint> evalFunc;
        jacobianData.evaluateReverse(start, end, evalFunc, tape, data);
      }

      /// Performs the AD \ref sec_forwardAD "forward" equation for a statement.
      template<typename Adjoint>
      CODI_INLINE static void incrementTangents(Adjoint const* const adjointVector, Adjoint& lhsAdjoint,
                                                Config::ArgumentSize const& numberOfArguments, size_t& curJacobianPos,
                                                Real const* const rhsJacobians,
                                                Identifier const* const rhsIdentifiers) {
        size_t endJacobianPos = curJacobianPos + numberOfArguments;

        while (curJacobianPos < endJacobianPos) CODI_Likely {
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
                                                                              NestedPosition const& end, Impl& tape,
                                                                              Adjoint* data,
                                                                              JacobianData& jacobianData) {
        Wrap_internalEvaluateForward_Step3_EvalStatements<Adjoint> evalFunc;
        jacobianData.evaluateForward(start, end, evalFunc, tape, data);
      }

    public:

      /// @name Functions from CustomAdjointVectorEvaluationTapeInterface
      /// @{

      using Base::evaluate;

      /// \copydoc codi::CustomAdjointVectorEvaluationTapeInterface::evaluate()
      template<typename Adjoint>
      CODI_NO_INLINE void evaluate(Position const& start, Position const& end, Adjoint* data) {
        VectorAccess<Adjoint> adjointWrapper(data);

        EventSystem<Impl>::notifyTapeEvaluateListeners(
            cast(), start, end, &adjointWrapper, EventHints::EvaluationKind::Reverse, EventHints::Endpoint::Begin);

        Base::internalEvaluateReverse_Step1_ExtFunc(
            start, end, JacobianBaseTape::template internalEvaluateReverse_Step2_DataExtraction<Adjoint>,
            &adjointWrapper, cast(), data, jacobianData);

        EventSystem<Impl>::notifyTapeEvaluateListeners(cast(), start, end, &adjointWrapper,
                                                       EventHints::EvaluationKind::Reverse, EventHints::Endpoint::End);
      }

      /// \copydoc codi::CustomAdjointVectorEvaluationTapeInterface::evaluate()
      template<typename Adjoint>
      CODI_NO_INLINE void evaluateForward(Position const& start, Position const& end, Adjoint* data) {
        VectorAccess<Adjoint> adjointWrapper(data);

        EventSystem<Impl>::notifyTapeEvaluateListeners(
            cast(), start, end, &adjointWrapper, EventHints::EvaluationKind::Forward, EventHints::Endpoint::Begin);

        Base::internalEvaluateForward_Step1_ExtFunc(
            start, end, JacobianBaseTape::template internalEvaluateForward_Step2_DataExtraction<Adjoint>,
            &adjointWrapper, cast(), data, jacobianData);

        EventSystem<Impl>::notifyTapeEvaluateListeners(cast(), start, end, &adjointWrapper,
                                                       EventHints::EvaluationKind::Forward, EventHints::Endpoint::End);
      }

      /// @}
      /*******************************************************************************/
      /// @name Functions from DataManagementTapeInterface
      /// @{

      /// \copydoc codi::DataManagementTapeInterface::swap()
      CODI_INLINE void swap(Impl& other) {
        // Index manager does not need to be swapped, it is either static or swapped with the vector data.
        // Vectors are swapped recursively in the base class.

        adjoints.swap(other.adjoints);

        Base::swap(other);
      }

      /// \copydoc codi::DataManagementTapeInterface::deleteAdjointVector()
      void deleteAdjointVector() {
        adjoints.resize(1);
      }

      /// \copydoc codi::DataManagementTapeInterface::resizeAdjointVector()
      void resizeAdjointVector() {
        checkAdjointSize(indexManager.get().getLargestCreatedIndex());
      }

      /// \copydoc codi::DataManagementTapeInterface::beginUseAdjointVector()
      void beginUseAdjointVector() {
        adjoints.beginUse();
      }

      /// \copydoc codi::DataManagementTapeInterface::endUseAdjointVector()
      void endUseAdjointVector() {
        adjoints.endUse();
      }

      /// \copydoc codi::DataManagementTapeInterface::getParameter()
      size_t getParameter(TapeParameters parameter) const {
        switch (parameter) {
          case TapeParameters::AdjointSize:
            return adjoints.size();
            break;
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
          case TapeParameters::AdjointSize:
            adjoints.resize(value);
            break;
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

      /// \copydoc codi::DataManagementTapeInterface::createVectorAccess()
      VectorAccess<Gradient>* createVectorAccess() {
        return createVectorAccessCustomAdjoints(adjoints.data());
      }

      /// \copydoc codi::DataManagementTapeInterface::createVectorAccessCustomAdjoints()
      template<typename Adjoint>
      VectorAccess<Adjoint>* createVectorAccessCustomAdjoints(Adjoint* data) {
        return new VectorAccess<Adjoint>(data);
      }

      /// \copydoc codi::DataManagementTapeInterface::deleteVectorAccess()
      void deleteVectorAccess(VectorAccessInterface<Real, Identifier>* access) {
        delete access;
      }

      /// @}
      /*******************************************************************************/
      /// @name Functions from ExternalFunctionTapeInterface
      /// @{

      /// \copydoc codi::ExternalFunctionTapeInterface::registerExternalFunctionOutput()
      template<typename Lhs>
      Real registerExternalFunctionOutput(LhsExpressionInterface<Real, Gradient, Impl, Lhs>& value) {
        cast().internalRegisterInput(value, false);

        return Real();
      }

      /// @}
      /*******************************************************************************/
      /// @name Functions from ForwardEvaluationTapeInterface
      /// @{

      using Base::evaluateForward;

      /// \copydoc codi::ForwardEvaluationTapeInterface::evaluateForward()
      void evaluateForward(Position const& start, Position const& end,
                           AdjointsManagement adjointsManagement = AdjointsManagement::Automatic) {
        if (AdjointsManagement::Automatic == adjointsManagement) {
          checkAdjointSize(indexManager.get().getLargestCreatedIndex());
          adjoints.beginUse();
        }

        codiAssert(indexManager.get().getLargestCreatedIndex() < (Identifier)adjoints.size());

        cast().evaluateForward(start, end, adjoints.data());

        if (AdjointsManagement::Automatic == adjointsManagement) {
          adjoints.endUse();
        }
      }

      /// @}
      /*******************************************************************************/
      /// @name Functions from ManualStatementPushTapeInterface
      /// @{

      /// \copydoc codi::ManualStatementPushTapeInterface::pushJacobianManual()
      void pushJacobianManual(Real const& jacobian, Real const& value, Identifier const& index) {
        CODI_UNUSED(value);

        cast().incrementManualPushCounter();

        jacobianData.pushData(jacobian, index);

        if (Config::StatementEvents) {
          if (this->manualPushCounter == this->manualPushGoal) {
            // emit statement event
            Real* jacobians;
            Identifier* rhsIdentifiers;
            jacobianData.getDataPointers(jacobianData.reserveItems(0), jacobians, rhsIdentifiers);
            jacobians -= this->manualPushGoal;
            rhsIdentifiers -= this->manualPushGoal;

            EventSystem<Impl>::notifyStatementStoreOnTapeListeners(cast(), this->manualPushLhsIdentifier,
                                                                   this->manualPushLhsValue, this->manualPushGoal,
                                                                   rhsIdentifiers, jacobians);
          }
        }
      }

      /// \copydoc codi::ManualStatementPushTapeInterface::storeManual()
      void storeManual(Real const& lhsValue, Identifier& lhsIndex, Config::ArgumentSize const& size) {
        CODI_UNUSED(lhsValue);

        codiAssert(size < Config::MaxArgumentSize);

        statementData.reserveItems(1);
        jacobianData.reserveItems(size);

        indexManager.get().template assignIndex<Impl>(lhsIndex);
        cast().pushStmtData(lhsIndex, (Config::ArgumentSize)size);

        cast().initializeManualPushData(lhsValue, lhsIndex, size);
      }

      /// @}
      /*******************************************************************************/
      /// @name Functions from PositionalEvaluationTapeInterface
      /// @{

      /// \copydoc codi::PositionalEvaluationTapeInterface::evaluate()
      CODI_INLINE void evaluate(Position const& start, Position const& end,
                                AdjointsManagement adjointsManagement = AdjointsManagement::Automatic) {
        if (AdjointsManagement::Automatic == adjointsManagement) {
          checkAdjointSize(indexManager.get().getLargestCreatedIndex());
          adjoints.beginUse();
        }

        codiAssert(indexManager.get().getLargestCreatedIndex() < (Identifier)adjoints.size());

        evaluate(start, end, adjoints.data());

        if (AdjointsManagement::Automatic == adjointsManagement) {
          adjoints.endUse();
        }
      }

      /// @}
      /*******************************************************************************/
      /// @name Functions from PreaccumulationEvaluationTapeInterface
      /// @{

      /// \copydoc codi::PreaccumulationEvaluationTapeInterface::evaluateKeepState()
      void evaluateKeepState(Position const& start, Position const& end,
                             AdjointsManagement adjointsManagement = AdjointsManagement::Automatic) {
        evaluate(start, end, adjointsManagement);
      }

      /// \copydoc codi::PreaccumulationEvaluationTapeInterface::evaluateForwardKeepState()
      void evaluateForwardKeepState(Position const& start, Position const& end,
                                    AdjointsManagement adjointsManagement = AdjointsManagement::Automatic) {
        evaluateForward(start, end, adjointsManagement);
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

    private:

      CODI_INLINE void checkAdjointSize(Identifier const& identifier) {
        if (identifier >= (Identifier)adjoints.size()) {
          internalResizeAdjointsVector();
        }
      }

      CODI_NO_INLINE void internalResizeAdjointsVector() {
        // overallocate as next multiple of Config::ChunkSize
        adjoints.resize(getNextMultiple((size_t)indexManager.get().getLargestCreatedIndex() + 1, Config::ChunkSize));
      }
  };
}
