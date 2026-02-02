/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2026 Chair for Scientific Computing (SciComp), University of Kaiserslautern-Landau
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

#include <algorithm>
#include <cmath>
#include <functional>
#include <type_traits>
#include <utility>

#include "../config.h"
#include "../expressions/aggregate/aggregatedActiveType.hpp"
#include "../expressions/aggregate/arrayAccessExpression.hpp"
#include "../expressions/emptyExpression.hpp"
#include "../expressions/lhsExpressionInterface.hpp"
#include "../expressions/logic/compileTimeTraversalLogic.hpp"
#include "../expressions/logic/constructStaticContext.hpp"
#include "../expressions/logic/helpers/forEachLeafLogic.hpp"
#include "../expressions/logic/helpers/jacobianComputationLogic.hpp"
#include "../expressions/logic/helpers/mathStatementGenLogic.hpp"
#include "../expressions/logic/traversalLogic.hpp"
#include "../misc/demangleName.hpp"
#include "../misc/macros.hpp"
#include "../misc/mathUtility.hpp"
#include "../misc/memberStore.hpp"
#include "../traits/expressionTraits.hpp"
#include "commonTapeImplementation.hpp"
#include "data/chunk.hpp"
#include "data/chunkedData.hpp"
#include "indices/indexManagerInterface.hpp"
#include "misc/assignStatement.hpp"
#include "misc/primalAdjointVectorAccess.hpp"
#include "statementEvaluators/statementEvaluatorInterface.hpp"
#include "statementEvaluators/statementEvaluatorTapeInterface.hpp"
/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Type definitions for the primal value tapes.
   *
   * @tparam T_Real                See TapeTypesInterface.
   * @tparam T_Gradient            See TapeTypesInterface.
   * @tparam T_IndexManager        Index manager for the tape. Needs to implement IndexManagerInterface.
   * @tparam T_StatementEvaluator  Statement handle generator. Needs to implement StatementEvaluatorInterface and
   *                              StatementEvaluatorInnerTapeInterface.
   * @tparam T_Data                See TapeTypesInterface.
   */
  template<typename T_Real, typename T_Gradient, typename T_IndexManager, typename T_StatementEvaluator,
           template<typename, typename> class T_Data>
  struct PrimalValueTapeTypes : public TapeTypesInterface {
    public:

      using Real = CODI_DD(T_Real, double);                                              ///< See PrimalValueTapeTypes.
      using Gradient = CODI_DD(T_Gradient, double);                                      ///< See PrimalValueTapeTypes.
      using IndexManager = CODI_DD(T_IndexManager, CODI_T(IndexManagerInterface<int>));  ///< See PrimalValueTapeTypes.
      using StatementEvaluator = CODI_DD(CODI_T(T_StatementEvaluator),
                                         CODI_T(StatementEvaluatorInterface));  ///< See PrimalValueTapeTypes.
      template<typename Chunk, typename Nested>
      using Data = CODI_DD(CODI_T(T_Data<Chunk, Nested>),
                           CODI_T(DataInterface<Nested>));  ///< See PrimalValueTapeTypes.

      using Identifier = typename IndexManager::Index;  ///< See IndexManagerInterface.
      using ActiveTypeTapeData =
          typename IndexManager::ActiveTypeIndexData;     ///< Take the active real data from the index manager.
      using PassiveReal = RealTraits::PassiveReal<Real>;  ///< Basic computation type.

      constexpr static bool IsLinearIndexHandler = IndexManager::IsLinear;  ///< True if the index manager is linear.
      constexpr static bool IsStaticIndexHandler =
          IndexManager::NeedsStaticStorage;  ///< True if the index manager must be stored statically in the tape.

      using EvalHandle = typename StatementEvaluator::Handle;  ///< Handle type returned by the statement generator.

      /// Statement chunk contains \<argument size, eval handle, data size of statement in bytes\>.
      using StatementChunk = Chunk3<Config::ArgumentSize, EvalHandle, Config::LowLevelFunctionDataSize>;
      using StatementData = Data<StatementChunk, IndexManager>;  ///< Statement data vector.

      using StatementByteChunk = Chunk1<char>;                            ///< Binary data for the statements.
      using StatementByteData = Data<StatementByteChunk, StatementData>;  ///< Statement byte data vector.

      using NestedData = StatementByteData;  ///< See TapeTypesInterface.
  };

  /**
   * @brief Base class for all standard Primal value tape implementations.
   *
   * This class provides nearly a full implementation of the FullTapeInterface. There are just a few internal methods
   * left which need to be implemented by the final classes. These methods depend significantly on the index management
   * scheme and are performance critical.
   *
   * Tape evaluations are performed in several steps. Each methods calls the next method:
   * - evaluate
   * - internalEvaluate*
   * - internalEvaluate*_EvalStatements
   * The placeholder stands for Reverse, Forward, or Primal.
   *
   * @tparam T_TapeTypes needs to implement PrimalValueTapeTypes.
   * @tparam T_Impl Type of the final implementations
   */
  template<typename T_TapeTypes, typename T_Impl>
  struct PrimalValueBaseTape : public CommonTapeImplementation<T_TapeTypes, T_Impl>,
                               public StatementEvaluatorTapeInterface,
                               public StatementEvaluatorInnerTapeInterface {
    public:

      /// See PrimalValueBaseTape.
      using TapeTypes = CODI_DD(T_TapeTypes,
                                CODI_T(PrimalValueTapeTypes<double, double, IndexManagerInterface<int>,
                                                            StatementEvaluatorInterface, DefaultChunkedData>));
      /// See PrimalValueBaseTape.
      using Impl = CODI_DD(T_Impl, CODI_T(PrimalValueBaseTape));

      using Base = CommonTapeImplementation<T_TapeTypes, T_Impl>;  ///< Base class abbreviation.
      friend Base;  ///< Allow the base class to call protected and private methods.

      using Real = typename TapeTypes::Real;                              ///< See TapeTypesInterface.
      using Gradient = typename TapeTypes::Gradient;                      ///< See TapeTypesInterface.
      using IndexManager = typename TapeTypes::IndexManager;              ///< See TapeTypesInterface.
      using StatementEvaluator = typename TapeTypes::StatementEvaluator;  ///< See PrimalValueTapeTypes.
      using Identifier = typename TapeTypes::Identifier;                  ///< See PrimalValueTapeTypes.
      using ActiveTypeTapeData = typename TapeTypes::ActiveTypeTapeData;  ///< See TapeTypesInterface.

      using EvalHandle = typename TapeTypes::EvalHandle;  ///< See PrimalValueTapeTypes.

      using StatementData = typename TapeTypes::StatementData;          ///< See PrimalValueTapeTypes.
      using StatementByteData = typename TapeTypes::StatementByteData;  ///< See PrimalValueTapeTypes.

      using PassiveReal = RealTraits::PassiveReal<Real>;  ///< Basic computation type.

      using NestedPosition = typename StatementData::Position;  ///< See PrimalValueTapeTypes.
      using Position = typename Base::Position;                 ///< See TapeTypesInterface.

      /// Vector access type generated by this tape.
      template<typename AdjointVector>
      using VectorAccess = PrimalAdjointVectorAccess<Real, Identifier, AdjointVector>;

      static bool constexpr AllowJacobianOptimization = false;  ///< See InternalStatementRecordingTapeInterface.
      static bool constexpr HasPrimalValues = true;             ///< See PrimalEvaluationTapeInterface.
      static bool constexpr LinearIndexHandling =
          TapeTypes::IsLinearIndexHandler;  ///< See IdentifierInformationTapeInterface.
      static bool constexpr RequiresPrimalRestore =
          !TapeTypes::IsLinearIndexHandler;  ///< See PrimalEvaluationTapeInterface.

      /// Used for generating arrays for lhs handling.
      template<typename T>
      using StackArray = std::array<T, Config::MaxArgumentSize>;

    protected:

      ///< Expressions for manual statement pushes.
      static EvalHandle const jacobianExpressions[Config::MaxArgumentSize];

      MemberStore<IndexManager, Impl, TapeTypes::IsStaticIndexHandler> indexManager;  ///< Index manager.
      StatementData statementData;          ///< Data stream for statement specific data.
      StatementByteData statementByteData;  ///< Data stream for statement byte data.

      std::vector<Gradient> adjoints;  ///< Evaluation vector for AD.
      std::vector<Real> primals;       ///< Current state of primal values in the program.
      std::vector<Real> primalsCopy;   ///< Copy of primal values for AD evaluations.

      Real* manualPushJacobians;          ///< Stores the pointer to the array for the Jacobian values of a manual
                                          ///< statement push.
      Identifier* manualPushIdentifiers;  ///< Stores the pointer to the array for the Jacobian values of a manual
                                          ///< statement push.

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
      static void internalEvaluateForward_EvalStatements(Args&&... args);

      /// Perform a primal evaluation of the tape. Arguments are from the recursive eval methods of the DataInterface.
      template<typename... Args>
      static void internalEvaluatePrimal_EvalStatements(Args&&... args);

      /// Perform a reverse evaluation of the tape. Arguments are from the recursive eval methods of the DataInterface.
      template<typename... Args>
      static void internalEvaluateReverse_EvalStatements(Args&&... args);

      /// Reset the primal values to the given position.
      void internalResetPrimalValues(Position const& pos);

      /// @}

    public:

      /// Constructor
      PrimalValueBaseTape()
          : Base(),
            indexManager(Config::MaxArgumentSize),  // Reserve first items for passive values.
            statementData(Config::ChunkSize),
            statementByteData(Config::ByteDataChunkSize),
            adjoints(1),  // Ensure that adjoint[0] exists, see its use in gradient() const.
            primals(0),
            primalsCopy(0) {
        checkPrimalSize(true);

        statementData.setNested(&indexManager.get());
        statementByteData.setNested(&statementData);

        Base::init(&statementByteData);

        Base::options.insert(TapeParameters::AdjointSize);
        Base::options.insert(TapeParameters::LargestIdentifier);
        Base::options.insert(TapeParameters::PrimalSize);
        Base::options.insert(TapeParameters::StatementSize);
        Base::options.insert(TapeParameters::StatementByteSize);
      }

      /*******************************************************************************/
      /// @name Functions from GradientAccessTapeInterface
      /// @{

      /// \copydoc codi::GradientAccessTapeInterface::gradient(T_Identifier const&, AdjointsManagement)
      /// <br> Implementation: Automatic adjoints management only involves bounds checking and resizing. Primal value
      /// tapes do not implement adjoints locking.
      CODI_INLINE Gradient& gradient(Identifier const& identifier,
                                     AdjointsManagement adjointsManagement = AdjointsManagement::Automatic) {
        if (AdjointsManagement::Automatic == adjointsManagement) {
          checkAdjointSize(identifier);
        }

        codiAssert(identifier < (Identifier)adjoints.size());

        return adjoints[identifier];
      }

      /// \copydoc codi::GradientAccessTapeInterface::gradient(T_Identifier const&, AdjointsManagement) const
      /// <br> Implementation: Automatic adjoints management only involves bounds checking. Primal value tapes do not
      /// implement adjoints locking.
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

      /// \copydoc codi::InternalStatementRecordingTapeInterface::initTapeData()
      template<typename Real>
      CODI_INLINE void initTapeData(Real& value, ActiveTypeTapeData& data) {
        CODI_UNUSED(value);

        indexManager.get().initIndex(data);
      }

      /// \copydoc codi::InternalStatementRecordingTapeInterface::destroyTapeData()
      template<typename Real>
      CODI_INLINE void destroyTapeData(Real& value, ActiveTypeTapeData& data) {
        CODI_UNUSED(value);

        indexManager.get().template freeIndex<Impl>(data);
      }

      /// @}

      /// Definition of the data for a statement. The pointers a populated from a byte data stream.
      struct StatementDataPointers {
          Identifier* CODI_RESTRICT rhsIdentifiers;   ///< Array for the rhs identifiers.
          Identifier* CODI_RESTRICT lhsIdentifiers;   ///< Array for the lhs identifiers.
          Real* CODI_RESTRICT passiveValues;          ///< Array for the passive values of a statement.
          Real* CODI_RESTRICT oldLhsValues;           ///< Array for the old values of the lhs values.
          PassiveReal* CODI_RESTRICT constantValues;  ///< Array for the constant values in the statement.

          /// Push data into an array. Sets the value at the current position and then changes the pointer to the next
          /// position.
          template<typename T>
          CODI_INLINE void push(T* CODI_RESTRICT& array, T const& value) {
            *array = value;
            array += 1;
          }

          /// Set all the pointer from the given byte data. The sizes are used to compute the offsets.
          CODI_INLINE void populate(size_t lhsSize, size_t rhsSize, size_t passiveSize, size_t constantSize,
                                    char* byteData) {
            size_t reserverLhsSize = LinearIndexHandling ? 0 : lhsSize;

            char* curPos = byteData;

            rhsIdentifiers = reinterpret_cast<Identifier*>(curPos);
            curPos += sizeof(Identifier) * rhsSize;

            lhsIdentifiers = reinterpret_cast<Identifier*>(curPos);
            curPos += sizeof(Identifier) * reserverLhsSize;

            passiveValues = reinterpret_cast<Real*>(curPos);
            curPos += sizeof(Real) * passiveSize;

            oldLhsValues = reinterpret_cast<Real*>(curPos);
            curPos += sizeof(Real) * reserverLhsSize;

            constantValues = reinterpret_cast<PassiveReal*>(curPos);
            curPos += sizeof(PassiveReal) * constantSize;

            codiAssert(curPos == byteData + computeSize(lhsSize, rhsSize, passiveSize, constantSize));
          }

          /// Computes the byte size for the data of the statement.
          CODI_INLINE Config::LowLevelFunctionDataSize computeSize(size_t lhsSize, size_t rhsSize, size_t passiveSize,
                                                                   size_t constantSize) {
            size_t reserverLhsSize = LinearIndexHandling ? 0 : lhsSize;

            return sizeof(Identifier) * (rhsSize + reserverLhsSize) + sizeof(Real) * (passiveSize + reserverLhsSize) +
                   sizeof(PassiveReal) * (constantSize);
          }
      };

    protected:

      /// Add data for a lhs entry to the data streams.
      CODI_INLINE void pushLhsData(Identifier const& index, Real const& oldPrimalValue,
                                   StatementDataPointers& pointers) {
        if (!LinearIndexHandling) {
          pointers.push(pointers.lhsIdentifiers, index);
          pointers.push(pointers.oldLhsValues, oldPrimalValue);
        }
      }

      /// Reserve all the data for a statement.
      CODI_INLINE Config::LowLevelFunctionDataSize reserveStmtDataManual(size_t lhsSize, size_t rhsSize,
                                                                         size_t passiveSize, size_t constantSize,
                                                                         StatementDataPointers& pointers) {
        Config::LowLevelFunctionDataSize byteSize = pointers.computeSize(lhsSize, rhsSize, passiveSize, constantSize);
        statementData.reserveItems(1);
        statementByteData.reserveItems(byteSize);

        char* byteData = nullptr;
        statementByteData.getDataPointers(byteData);
        statementByteData.addDataSize(byteSize);

        pointers.populate(lhsSize, rhsSize, passiveSize, constantSize, byteData);

        return byteSize;
      }

      /// Reserve all the data for a statement.
      template<typename Lhs, typename Rhs>
      CODI_INLINE Config::LowLevelFunctionDataSize reserveStmtData(size_t activeArguments,
                                                                   StatementDataPointers& pointers) {
        size_t constexpr MaxActiveArgs = ExpressionTraits::NumberOfActiveTypeArguments<Rhs>::value;
        size_t constexpr MaxConstantArgs = ExpressionTraits::NumberOfConstantTypeArguments<Rhs>::value;
        size_t constexpr MaxOutputArgs = ExpressionTraits::NumberOfActiveTypeArguments<Lhs>::value;

        return reserveStmtDataManual(MaxOutputArgs, MaxActiveArgs, MaxActiveArgs - activeArguments, MaxConstantArgs,
                                     pointers);
      }

      /// Count all arguments that have non-zero index.
      struct CountActiveArguments : public ForEachLeafLogic<CountActiveArguments> {
        public:

          /// \copydoc codi::ForEachLeafLogic::handleActive
          template<typename Node>
          CODI_INLINE void handleActive(Node const& node, size_t& numberOfActiveArguments, IndexManager& indexManager) {
            indexManager.validateRhsIndex(node.getTapeData());

            if (CODI_ENABLE_CHECK(Config::CheckZeroIndex, IndexManager::InactiveIndex != node.getIdentifier())) {
              numberOfActiveArguments += 1;
            }
          }
      };

      /// Push all data for each argument.
      struct PushIdentifierPassiveAndConstant : public ForEachLeafLogic<PushIdentifierPassiveAndConstant> {
        public:

          /// \copydoc codi::ForEachLeafLogic::handleActive
          template<typename Node>
          CODI_INLINE void handleActive(Node const& node, StatementDataPointers& pointers, size_t& curPassiveArgument) {
            Identifier rhsIndex = node.getIdentifier();
            if (CODI_ENABLE_CHECK(Config::CheckZeroIndex, IndexManager::InactiveIndex == rhsIndex)) {
              rhsIndex = curPassiveArgument;

              curPassiveArgument += 1;
              pointers.push(pointers.passiveValues, node.getValue());
            }

            pointers.push(pointers.rhsIdentifiers, rhsIndex);
          }

          /// \copydoc codi::ForEachLeafLogic::handleConstant
          template<typename Node>
          CODI_INLINE void handleConstant(Node const& node, StatementDataPointers& pointers,
                                          size_t& curPassiveArgument) {
            CODI_UNUSED(curPassiveArgument);

            using AggregatedTraits = codi::RealTraits::AggregatedTypeTraits<typename Node::Real>;
            using ConversionOperator = typename Node::template ConversionOperator<PassiveReal>;

            typename Node::Real v = node.getValue();

            static_for<AggregatedTraits::Elements>([&](auto i) CODI_LAMBDA_INLINE {
              pointers.push(pointers.constantValues,
                            ConversionOperator::toDataStore(AggregatedTraits::template arrayAccess<i.value>(v)));
            });
          }
      };

      /// Store all data for the rhs and the statement data. Does not do anything if there are no active arguments.
      /// Does not store the lhs data.
      ///
      /// @return true if the rhs data was stored and the lhs data needs to be stored.
      template<typename Lhs, typename RhsType, typename Rhs>
      CODI_INLINE bool storeArgumentAndStmtData(ExpressionInterface<RhsType, Rhs> const& rhs,
                                                StatementDataPointers& pointers) {
        using Stmt = AssignStatement<Lhs, Rhs>;

        CountActiveArguments countActiveArguments;
        PushIdentifierPassiveAndConstant pushStatement;

        codiAssert(ExpressionTraits::NumberOfActiveTypeArguments<Rhs>::value < Config::MaxArgumentSize);
        codiAssert(ExpressionTraits::NumberOfConstantTypeArguments<Rhs>::value < Config::MaxArgumentSize);
        codiAssert(ExpressionTraits::NumberOfActiveTypeArguments<Lhs>::value < Config::MaxArgumentSize);

        size_t activeArguments = 0;
        countActiveArguments.eval(rhs.cast(), activeArguments, indexManager.get());

        if (CODI_ENABLE_CHECK(Config::CheckEmptyStatements, 0 != activeArguments)) {
          Config::LowLevelFunctionDataSize byteSize = reserveStmtData<Lhs, Rhs>(activeArguments, pointers);

          size_t passiveArguments = 0;
          pushStatement.eval(rhs.cast(), pointers, passiveArguments);
          statementData.pushData((Config::ArgumentSize)passiveArguments,
                                 StatementEvaluator::template createHandle<Impl, Impl, Stmt>(), byteSize);

          return true;
        } else {
          return false;
        }
      }

      /// Computes Jacobian entries for the event system.
      struct JacobianExtractionLogic : public JacobianComputationLogic<JacobianExtractionLogic> {
        private:
          size_t pos;

        public:

          /// Constructor
          JacobianExtractionLogic() : pos(0) {}

          /// Stores the identifiers and Jacobians.
          template<typename Node, typename Jacobian>
          CODI_INLINE void handleJacobianOnActive(Node const& node, Jacobian jacobianExpr, Identifier* rhsIdentifiers,
                                                  Real* jacobians) {
            rhsIdentifiers[pos] = node.getIdentifier();
            jacobians[pos] = jacobianExpr;
            pos++;
          }
      };

    public:

      /// @{
      ///
      ///       /// \copydoc codi::InternalStatementRecordingTapeInterface::store() <br>
      /// Implementation for AggregatedActiveType.
      template<typename Aggregated, typename Type, typename Lhs, typename Rhs>
      CODI_INLINE void store(AggregatedActiveType<Aggregated, Type, Lhs>& lhs,
                             ExpressionInterface<Aggregated, Rhs> const& rhs) {
        using AggregatedTraits = RealTraits::AggregatedTypeTraits<Aggregated>;
        int constexpr Elements = AggregatedTraits::Elements;

        bool primalStored = false;
        StatementDataPointers pointers = {};

        if (storeArgumentAndStmtData<Lhs>(rhs, pointers)) {
          bool generatedNewIndex = false;
          static_for<Elements>([&](auto i) CODI_LAMBDA_INLINE {
            generatedNewIndex |= indexManager.get().template assignIndex<Impl>(lhs.values[i.value].getTapeData());
          });
          checkPrimalSize(generatedNewIndex);

          Aggregated real = rhs.cast().getValue();
          static_for<Elements>([&](auto i) CODI_LAMBDA_INLINE {
            Identifier lhsIdentifier = lhs.values[i.value].getIdentifier();
            Real& primalEntry = primals[lhsIdentifier];

            pushLhsData(lhsIdentifier, primalEntry, pointers);

            if (Config::StatementEvents) {
              size_t constexpr MaxActiveArgs = ExpressionTraits::NumberOfActiveTypeArguments<Rhs>::value;

              JacobianExtractionLogic getRhsIdentifiersAndJacobians;
              std::array<Identifier, MaxActiveArgs> rhsIdentifiers;
              std::array<Real, MaxActiveArgs> jacobians;
              getRhsIdentifiersAndJacobians.eval(ArrayAccessExpression<Aggregated, i.value, Rhs>(rhs), Real(1.0),
                                                 rhsIdentifiers.data(), jacobians.data());

              EventSystem<Impl>::notifyStatementStoreOnTapeListeners(
                  cast(), lhsIdentifier, AggregatedTraits::template arrayAccess<i.value>(real), MaxActiveArgs,
                  rhsIdentifiers.data(), jacobians.data());
            }

            lhs.values[i.value].value() = AggregatedTraits::template arrayAccess<i.value>(real);
            primalEntry = lhs.values[i.value].getValue();
          });

          primalStored = true;
        }

        if (!primalStored) {
          Aggregated real = rhs.cast().getValue();

          static_for<Elements>([&](auto i) CODI_LAMBDA_INLINE {
            lhs.values[i.value].value() = AggregatedTraits::template arrayAccess<i.value>(real);
            indexManager.get().template freeIndex<Impl>(lhs.values[i.value].getTapeData());
          });
        }
      }

      /// \copydoc codi::InternalStatementRecordingTapeInterface::store() <br>
      /// Optimization for copy statements of aggregated types.
      template<typename Aggregated, typename Type, typename Lhs, typename Rhs>
      CODI_INLINE void store(AggregatedActiveType<Aggregated, Type, Lhs>& lhs,
                             AggregatedActiveType<Aggregated, Type, Rhs> const& rhs) {
        using AggregatedTraits = RealTraits::AggregatedTypeTraits<Aggregated>;

        int constexpr Elements = AggregatedTraits::Elements;

        if (CODI_ENABLE_CHECK(Config::CheckTapeActivity, cast().isActive())) {
          if (IndexManager::CopyNeedsStatement || !Config::CopyOptimization) {
            store(lhs, static_cast<ExpressionInterface<Aggregated, Rhs> const&>(rhs));
            return;
          } else {
            static_for<Elements>([&](auto i) CODI_LAMBDA_INLINE {
              indexManager.get().template copyIndex<Impl>(lhs.values[i.value].getTapeData(),
                                                          rhs.values[i.value].getTapeData());
            });
          }
        } else {
          static_for<Elements>([&](auto i) CODI_LAMBDA_INLINE {
            indexManager.get().template freeIndex<Impl>(lhs.values[i.value].getTapeData());
          });
        }

        static_for<Elements>([&](auto i) CODI_LAMBDA_INLINE {
          lhs.values[i.value].value() = rhs.values[i.value].getValue();
        });
      }

      /// \copydoc codi::InternalStatementRecordingTapeInterface::store()
      template<typename Lhs, typename Rhs>
      CODI_INLINE void store(LhsExpressionInterface<Real, Gradient, Impl, Lhs>& lhs,
                             ExpressionInterface<Real, Rhs> const& rhs) {
        bool primalStored = false;
        StatementDataPointers pointers = {};

        if (storeArgumentAndStmtData<Lhs>(rhs, pointers)) {
          bool generatedNewIndex = indexManager.get().template assignIndex<Impl>(lhs.cast().getTapeData());
          checkPrimalSize(generatedNewIndex);

          Real& primalEntry = primals[lhs.cast().getIdentifier()];
          pushLhsData(lhs.cast().getIdentifier(), primalEntry, pointers);

          primalEntry = rhs.cast().getValue();

          if (Config::StatementEvents) {
            size_t constexpr MaxActiveArgs = ExpressionTraits::NumberOfActiveTypeArguments<Rhs>::value;

            JacobianExtractionLogic getRhsIdentifiersAndJacobians;
            std::array<Identifier, MaxActiveArgs> rhsIdentifiers;
            std::array<Real, MaxActiveArgs> jacobians;
            getRhsIdentifiersAndJacobians.eval(rhs.cast(), Real(1.0), rhsIdentifiers.data(), jacobians.data());

            EventSystem<Impl>::notifyStatementStoreOnTapeListeners(cast(), lhs.cast().getIdentifier(),
                                                                   rhs.cast().getValue(), MaxActiveArgs,
                                                                   rhsIdentifiers.data(), jacobians.data());
          }

          primalStored = true;
        }

        if (!primalStored) {
          indexManager.get().template freeIndex<Impl>(lhs.cast().getTapeData());
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
            indexManager.get().template copyIndex<Impl>(lhs.cast().getTapeData(), rhs.cast().getTapeData());
          }
        } else {
          indexManager.get().template freeIndex<Impl>(lhs.cast().getTapeData());
        }

        lhs.cast().value() = rhs.cast().getValue();
      }

      /// \copydoc codi::InternalStatementRecordingTapeInterface::store() <br>
      /// Specialization for passive assignments.
      template<typename Lhs>
      CODI_INLINE void store(LhsExpressionInterface<Real, Gradient, Impl, Lhs>& lhs, Real const& rhs) {
        indexManager.get().template freeIndex<Impl>(lhs.cast().getTapeData());

        lhs.cast().value() = rhs;
      }

      /// @}
      /*******************************************************************************
       * Protected helper function for ReverseTapeInterface
       */

    protected:

      /// Add a new input to the tape and update the primal value vector.
      template<typename Lhs>
      CODI_INLINE Real internalRegisterInput(LhsExpressionInterface<Real, Gradient, Impl, Lhs>& value,
                                             bool unusedIndex) {
        if (TapeTypes::IsLinearIndexHandler) {
          statementData.reserveItems(1);
        }

        bool generatedNewIndex;
        if (unusedIndex) {
          generatedNewIndex = indexManager.get().template assignUnusedIndex<Impl>(value.cast().getTapeData());
        } else {
          generatedNewIndex = indexManager.get().template assignIndex<Impl>(value.cast().getTapeData());
        }
        checkPrimalSize(generatedNewIndex);

        Real& primalEntry = primals[value.cast().getIdentifier()];
        if (TapeTypes::IsLinearIndexHandler) {
          statementData.pushData(
              Config::StatementInputTag,
              StatementEvaluator::template createHandle<Impl, Impl, AssignStatement<Lhs, EmptyExpression<Real>>>(), 0);
        }

        Real oldValue = primalEntry;
        primalEntry = value.cast().value();

        return oldValue;
      }

    public:

      /// @name Functions from ReverseTapeInterface
      /// @{

      /// \copydoc codi::ReverseTapeInterface::registerInput()
      template<typename Lhs>
      CODI_INLINE void registerInput(LhsExpressionInterface<Real, Gradient, Impl, Lhs>& value) {
        internalRegisterInput(value, true);
        EventSystem<Impl>::notifyTapeRegisterInputListeners(cast(), value.cast().value(), value.cast().getIdentifier());
      }

      /// \copydoc codi::ReverseTapeInterface::clearAdjoints(AdjointsManagement)
      /// <br> Implementation: Automatic adjoints management has no effect. Primal value tapes do not implement adjoints
      /// locking.
      CODI_INLINE void clearAdjoints(AdjointsManagement adjointsManagement = AdjointsManagement::Automatic) {
        CODI_UNUSED(adjointsManagement);

        size_t maxSize = std::min((size_t)indexManager.get().getLargestCreatedIndex() + 1, adjoints.size());
        for (size_t i = 0; i < maxSize; i += 1) {
          adjoints[i] = Gradient();
        }
      }

      /// \copydoc codi::ReverseTapeInterface::reset()
      /// <br> Implementation: Automatic adjoints management has no effect. Primal value tapes do not implement adjoints
      /// locking.
      CODI_INLINE void reset(bool resetAdjoints = true,
                             AdjointsManagement adjointsManagement = AdjointsManagement::Automatic) {
        for (Real& primal : primals) {
          primal = Real();
        }

        Base::reset(resetAdjoints, adjointsManagement);
      }

      /// @}

    protected:

      /// Adds data from all streams, the size of the adjoint vector, the size of the primal vector, and index manager
      /// data.
      CODI_INLINE TapeValues internalGetTapeValues() const {
        std::string name;
        if (TapeTypes::IsLinearIndexHandler) {
          name = "CoDi Tape Statistics ( PrimalValueLinearTape )";
        } else {
          name = "CoDi Tape Statistics ( PrimalValueReuseTape )";
        }
        TapeValues values = TapeValues(name);

        size_t nAdjoints = indexManager.get().getLargestCreatedIndex();
        double memoryAdjoints = static_cast<double>(nAdjoints) * static_cast<double>(sizeof(Gradient));

        size_t nPrimals = indexManager.get().getLargestCreatedIndex();
        double memoryPrimals = static_cast<double>(nPrimals) * static_cast<double>(sizeof(Real));

        values.addSection("Adjoint vector");
        values.addUnsignedLongEntry("Number of adjoints", nAdjoints);
        values.addDoubleEntry("Memory allocated", memoryAdjoints, TapeValues::LocalReductionOperation::Sum, true, true);

        values.addSection("Primal vector");
        values.addUnsignedLongEntry("Number of primals", nPrimals);
        values.addDoubleEntry("Memory allocated", memoryPrimals, TapeValues::LocalReductionOperation::Sum, true, true);

        values.addSection("Index manager");
        indexManager.get().addToTapeValues(values);

        values.addSection("Statement entries");
        statementData.addToTapeValues(values);

        values.addSection("Statement byte entries");
        statementByteData.addToTapeValues(values);

        return values;
      }

      /******************************************************************************
       * Protected helper function for CustomAdjointVectorEvaluationTapeInterface
       */

      /// Select the configured adjoint vector, see codi::Config::VariableAdjointInterfaceInPrimalTapes.
      template<typename AdjointVector>
      ADJOINT_VECTOR_TYPE* selectAdjointVector(VectorAccess<AdjointVector>* vectorAccess, AdjointVector data) {
        CODI_UNUSED(vectorAccess, data);

#if CODI_VariableAdjointInterfaceInPrimalTapes
        return vectorAccess;
#else
        CODI_STATIC_ASSERT(CODI_T(std::is_same<typename std::remove_reference<AdjointVector>::type, Gradient*>::value),
                           "Please enable 'CODI_VariableAdjointInterfaceInPrimalTapes' in order"
                           " to use custom adjoint vectors in the primal value tapes.");

        return data;
#endif
      }

      /// Perform the adjoint update based on the configuration in codi::Config::VariableAdjointInterfaceInPrimalTapes.
      struct IncrementForwardLogic : public JacobianComputationLogic<IncrementForwardLogic> {
        public:

          /// \copydoc codi::JacobianComputationLogic::handleJacobianOnActive()
          template<typename Node>
          CODI_INLINE void handleJacobianOnActive(Node const& node, Real jacobian, Gradient& lhsTangent,
                                                  ADJOINT_VECTOR_TYPE* adjointVector) {
            CODI_UNUSED(lhsTangent);

            if (CODI_ENABLE_CHECK(Config::IgnoreInvalidJacobians, RealTraits::isTotalFinite(jacobian))) CODI_Likely {
#if CODI_VariableAdjointInterfaceInPrimalTapes
              adjointVector->updateTangentWithLhs(node.getIdentifier(), jacobian);
#else
              lhsTangent += jacobian * adjointVector[node.getIdentifier()];
#endif
            }
          }
      };

      /// Additional wrapper that triggers compiler optimizations.
      CODI_WRAP_FUNCTION(Wrap_internalEvaluateForward_EvalStatements, Impl::internalEvaluateForward_EvalStatements);

      /// Internal method for the forward evaluation of the whole tape.
      template<bool copyPrimal, typename AdjointVector>
      CODI_NO_INLINE void internalEvaluateForward(Position const& start, Position const& end, AdjointVector&& data) {
        CODI_STATIC_ASSERT(
            Config::VariableAdjointInterfaceInPrimalTapes ||
                CODI_T(std::is_same<typename std::remove_reference<AdjointVector>::type, Gradient*>::value),
            "Please enable 'CODI_VariableAdjointInterfaceInPrimalTapes' in order"
            " to use custom adjoint vectors in the primal value tapes.");

        std::vector<Real> primalsCopy(0);
        Real* primalData = primals.data();

        if (copyPrimal) {
          primalsCopy = primals;
          primalData = primalsCopy.data();
        }

        VectorAccess<AdjointVector> vectorAccess(data, primalData);

        ADJOINT_VECTOR_TYPE* dataVector = selectAdjointVector<AdjointVector>(&vectorAccess, data);

        EventSystem<Impl>::notifyTapeEvaluateListeners(
            cast(), start, end, &vectorAccess, EventHints::EvaluationKind::Forward, EventHints::Endpoint::Begin);

        Wrap_internalEvaluateForward_EvalStatements evalFunc{};
        Base::llfByteData.evaluateForward(start, end, evalFunc, cast(), primalData, dataVector);

        EventSystem<Impl>::notifyTapeEvaluateListeners(cast(), start, end, &vectorAccess,
                                                       EventHints::EvaluationKind::Forward, EventHints::Endpoint::End);
      }

      /// Perform the adjoint update based on the configuration in codi::Config::VariableAdjointInterfaceInPrimalTapes.
      struct IncrementReversalLogic : public JacobianComputationLogic<IncrementReversalLogic> {
        public:

          /// See IncrementReversalLogic.
          template<typename Node>
          CODI_INLINE void handleJacobianOnActive(Node const& node, Real jacobian, Gradient const& lhsAdjoint,
                                                  ADJOINT_VECTOR_TYPE* adjointVector) {
            CODI_UNUSED(lhsAdjoint);

            if (CODI_ENABLE_CHECK(Config::IgnoreInvalidJacobians, RealTraits::isTotalFinite(jacobian))) CODI_Likely {
#if CODI_VariableAdjointInterfaceInPrimalTapes
              adjointVector->updateAdjointWithLhs(node.getIdentifier(), jacobian);
#else
              adjointVector[node.getIdentifier()] += jacobian * lhsAdjoint;
#endif
            }
          }
      };

      /// Additional wrapper that triggers compiler optimizations.
      CODI_WRAP_FUNCTION(Wrap_internalEvaluateReverse_EvalStatements, Impl::internalEvaluateReverse_EvalStatements);

      /// Internal method for the reverse evaluation of the whole tape.
      template<bool copyPrimal, typename AdjointVector>
      CODI_INLINE void internalEvaluateReverse(Position const& start, Position const& end, AdjointVector&& data) {
        CODI_STATIC_ASSERT(
            Config::VariableAdjointInterfaceInPrimalTapes ||
                CODI_T(std::is_same<typename std::remove_reference<AdjointVector>::type, Gradient*>::value),
            "Please enable 'CODI_VariableAdjointInterfaceInPrimalTapes' in order"
            " to use custom adjoint vectors in the primal value tapes.");

        Real* primalData = primals.data();

        if (copyPrimal) {
          primalsCopy = primals;
          primalData = primalsCopy.data();
        }

        VectorAccess<AdjointVector> vectorAccess(data, primalData);

        ADJOINT_VECTOR_TYPE* dataVector = selectAdjointVector<AdjointVector>(&vectorAccess, data);

        EventSystem<Impl>::notifyTapeEvaluateListeners(
            cast(), start, end, &vectorAccess, EventHints::EvaluationKind::Reverse, EventHints::Endpoint::Begin);

        Wrap_internalEvaluateReverse_EvalStatements evalFunc;
        Base::llfByteData.evaluateReverse(start, end, evalFunc, cast(), primalData, dataVector);

        EventSystem<Impl>::notifyTapeEvaluateListeners(cast(), start, end, &vectorAccess,
                                                       EventHints::EvaluationKind::Reverse, EventHints::Endpoint::End);
      }

    public:

      /// @name Functions from CustomAdjointVectorEvaluationTapeInterface
      /// @{

      using Base::evaluate;

      /// \copydoc codi::CustomAdjointVectorEvaluationTapeInterface::evaluate()
      template<typename AdjointVector>
      CODI_INLINE void evaluate(Position const& start, Position const& end, AdjointVector&& data) {
        internalEvaluateReverse<!TapeTypes::IsLinearIndexHandler>(start, end, std::forward<AdjointVector>(data));
      }

      /// \copydoc codi::CustomAdjointVectorEvaluationTapeInterface::evaluateForward()
      template<typename AdjointVector>
      CODI_INLINE void evaluateForward(Position const& start, Position const& end, AdjointVector&& data) {
        internalEvaluateForward<!TapeTypes::IsLinearIndexHandler>(start, end, std::forward<AdjointVector>(data));
      }

      /// \copydoc codi::CustomAdjointVectorEvaluationTapeInterface::getInternalAdjoints()
      CODI_INLINE Gradient* getInternalAdjoints() {
        return adjoints.data();
      }

      /// @}
      /*******************************************************************************/
      /// @name Functions from DataManagementTapeInterface
      /// @{

      /// \copydoc codi::DataManagementTapeInterface::swap()
      CODI_INLINE void swap(Impl& other) {
        // Index manager does not need to be swapped, it is either static or swapped with the vector data.
        // Vectors are swapped recursively in the base class.

        std::swap(adjoints, other.adjoints);
        std::swap(primals, other.primals);

        Base::swap(other);

        // Ensure that the primals vector of both tapes are sized according to the index manager.
        checkPrimalSize(true);
        other.checkPrimalSize(true);
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
      /// <br> Implementation: Empty since primal value tapes do not implement adjoints locking.
      void beginUseAdjointVector() {}

      /// \copydoc codi::DataManagementTapeInterface::endUseAdjointVector()
      /// <br> Implementation: Empty since primal value tapes do not implement adjoints locking.
      void endUseAdjointVector() {}

      /// \copydoc codi::DataManagementTapeInterface::getParameter()
      size_t getParameter(TapeParameters parameter) const {
        switch (parameter) {
          case TapeParameters::AdjointSize:
            return adjoints.size();
            break;
          case TapeParameters::LargestIdentifier:
            return indexManager.get().getLargestCreatedIndex();
            break;
          case TapeParameters::PrimalSize:
            return primals.size();
            break;
          case TapeParameters::StatementSize:
            return statementData.getDataSize();
            break;
          case TapeParameters::StatementByteSize:
            return statementByteData.getDataSize();
            break;
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
          case TapeParameters::LargestIdentifier:
            CODI_EXCEPTION("Tried to set a get only option.");
            break;
          case TapeParameters::PrimalSize:
            primals.resize(value);
            break;
          case TapeParameters::StatementSize:
            return statementData.resize(value);
            break;
          case TapeParameters::StatementByteSize:
            return statementByteData.resize(value);
            break;
          default:
            Base::setParameter(parameter, value);
            break;
        }
      }

      /// \copydoc codi::DataManagementTapeInterface::createVectorAccess()
      VectorAccess<Gradient*>* createVectorAccess() {
        return createVectorAccessCustomAdjoints(adjoints.data());
      }

      /// \copydoc codi::DataManagementTapeInterface::createVectorAccessCustomAdjoints()
      template<typename AdjointVector>
      VectorAccess<AdjointVector>* createVectorAccessCustomAdjoints(AdjointVector&& data) {
        return new VectorAccess<AdjointVector>(data, primals.data());
      }

      /// \copydoc codi::DataManagementTapeInterface::createVectorAccessCustomAdjoints()
      /// <br> Overload for pointers passed as lvalues. Ensures that the pointer is copied, not referenced.
      template<typename Adjoint>
      VectorAccess<Adjoint*>* createVectorAccessCustomAdjoints(Adjoint* data) {
        return new VectorAccess<Adjoint*>(data, primals.data());
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
        return internalRegisterInput(value, false);
      }

      /// @}
      /*******************************************************************************/
      /// @name Functions from ForwardEvaluationTapeInterface
      /// @{

      using Base::evaluateForward;

      /// \copydoc codi::ForwardEvaluationTapeInterface::evaluateForward()
      /// <br> Implementation: Automatic adjoints management only involves bounds checking and resizing. Primal value
      /// tapes do not implement adjoints locking.
      void evaluateForward(Position const& start, Position const& end,
                           AdjointsManagement adjointsManagement = AdjointsManagement::Automatic) {
        if (AdjointsManagement::Automatic == adjointsManagement) {
          checkAdjointSize(indexManager.get().getLargestCreatedIndex());
        }

        codiAssert(indexManager.get().getLargestCreatedIndex() < (Identifier)adjoints.size());

        cast().evaluateForward(start, end, adjoints.data());
      }

      /// @}
      /*******************************************************************************/
      /// @name Functions from ManualStatementPushTapeInterface
      /// @{

      /// \copydoc codi::ManualStatementPushTapeInterface::pushJacobianManual()
      void pushJacobianManual(Real const& jacobian, Real const& value, ActiveTypeTapeData const& data) {
        CODI_UNUSED(value);

        cast().incrementManualPushCounter();

        *manualPushIdentifiers = indexManager.get().getIndex(data);
        *manualPushJacobians = jacobian;

        manualPushIdentifiers += 1;
        manualPushJacobians += 1;

        if (Config::StatementEvents) {
          if (this->manualPushCounter == this->manualPushGoal) {
            // emit statement event
            manualPushIdentifiers -= this->manualPushGoal;
            manualPushJacobians -= this->manualPushGoal;

            EventSystem<Impl>::notifyStatementStoreOnTapeListeners(cast(), this->manualPushLhsIdentifier,
                                                                   this->manualPushLhsValue, this->manualPushGoal,
                                                                   manualPushIdentifiers, manualPushJacobians);
          }
        }
      }

      /// \copydoc codi::ManualStatementPushTapeInterface::storeManual()
      void storeManual(Real const& lhsValue, ActiveTypeTapeData& lhsData, Config::ArgumentSize const& size) {
        CODI_UNUSED(lhsValue);

        codiAssert(size < Config::MaxArgumentSize);

        StatementDataPointers pointers = {};
        Config::LowLevelFunctionDataSize byteSize = reserveStmtDataManual(1, size, size, 0, pointers);

        manualPushJacobians = pointers.passiveValues;
        manualPushIdentifiers = pointers.rhsIdentifiers;

        indexManager.get().template assignIndex<Impl>(lhsData);
        Real& primalEntry = primals[indexManager.get().getIndex(lhsData)];
        statementData.pushData(size, PrimalValueBaseTape::jacobianExpressions[size], byteSize);
        pushLhsData(indexManager.get().getIndex(lhsData), primalEntry, pointers);

        primalEntry = lhsValue;

        cast().initializeManualPushData(lhsValue, indexManager.get().getIndex(lhsData), size);
      }

      /// @}
      /*******************************************************************************/
      /// @name Functions from ReadWriteTapeInterface
      /// @{

      /// \copydoc codi::ReadWriteTapeInterface::createStatementManual(codi::ReadWriteTapeInterface::Identifier const&,
      /// codi::ReadWriteTapeInterface::Real const&, codi::Config::ArgumentSize const&,
      /// codi::ReadWriteTapeInterface::Identifier const* const, codi::Config::ArgumentSize const&,
      /// codi::ReadWriteTapeInterface::Real const* const, codi::Config::ArgumentSize const&,
      /// codi::ReadWriteTapeInterface::Real const* const, codi::ReadWriteTapeInterface::EvalHandle const&)
      void createStatementManual(Config::ArgumentSize const& nOutputValues, Identifier const* const lhsIndices,
                                 Real const* const lhsValues, Config::ArgumentSize const& nActiveValues,
                                 Identifier const* const rhsIdentifiers, Config::ArgumentSize const& nPassiveValues,
                                 Real const* const rhsPrimals, Config::ArgumentSize const& nConstants,
                                 Real const* const rhsConstant, EvalHandle const& evalHandle) {
        Impl& impl = cast();
        if (Config::StatementLowLevelFunctionTag == nPassiveValues) CODI_Unlikely {
          // TODO.
        } else if (Config::StatementInputTag == nPassiveValues && TapeTypes::IsLinearIndexHandler) CODI_Unlikely {
          statementData.reserveItems(1);
          primals[lhsIndices[0]] = lhsValues[0];
          statementData.pushData(Config::StatementInputTag, evalHandle, 0);
        } else CODI_Likely {
          StatementDataPointers pointers = {};
          Config::LowLevelFunctionDataSize byteSize =
              reserveStmtDataManual(nOutputValues, nActiveValues, nPassiveValues, nConstants, pointers);

          impl.initializeManualPushData(lhsValues[0], lhsIndices[0], nPassiveValues);

          for (size_t activeCount = 0; activeCount < nActiveValues; activeCount++) {
            pointers.rhsIdentifiers[activeCount] = rhsIdentifiers[activeCount];
          }
          for (size_t passiveCount = 0; passiveCount < nPassiveValues; passiveCount++) {
            cast().incrementManualPushCounter();
            pointers.passiveValues[passiveCount] = rhsPrimals[passiveCount];
          }
          for (size_t constantCount = 0; constantCount < nConstants; constantCount++) {
            pointers.constantValues[constantCount] = rhsConstant[constantCount];
          }

          statementData.pushData(nPassiveValues, evalHandle, byteSize);

          for (Config::ArgumentSize i = 0; i < nOutputValues; i += 1) {
            pushLhsData(lhsIndices[i], lhsValues[i], pointers);

            if (TapeTypes::IsLinearIndexHandler) {
              primals[lhsIndices[i]] = lhsValues[i];
            }
          }
        }
      }

      using Base::writeTape;

      /**
       * \copydoc codi::ReadWriteTapeInterface::writeTape(codi::ReadWriteTapeInterface::WriterInterface*,
       * codi:ReadWriteTapeInterface::Position const&, codi:ReadWriteTapeInterface::Position const&)
       */
      template<typename Type>
      CODI_NO_INLINE void writeTape(codi::TapeWriterInterface<Type>* writer, Position const& start,
                                    Position const& end) {
        Impl& impl = cast();
        writer->start(impl);
        Base::llfByteData.evaluateForward(start, end, Impl::template internalWriteTape<Type>, primals.data(), writer);
        writer->finish();
      }

      /// @}
      /*******************************************************************************/
      /// @name Functions from IdentifierInformationTapeInterface
      /// @{

      /// \copydoc codi::IdentifierInformationTapeInterface::getIndexManager()
      IndexManager& getIndexManager() {
        return indexManager.get();
      }

      /// @}
      /*******************************************************************************/
      /// @name Functions from LowLevelFunctionTapeInterface
      /// @{

      /// @copydoc LowLevelFunctionTapeInterface::pushLowLevelFunction
      CODI_INLINE void pushLowLevelFunction(Config::LowLevelFunctionToken token, size_t size, ByteDataView& data) {
        statementData.reserveItems(1);

        Base::internalStoreLowLevelFunction(token, size, data);

        statementData.pushData(Config::StatementLowLevelFunctionTag, EvalHandle(), 0);
      }

      /// @}
      /*******************************************************************************/
      /// @name Functions from PositionalEvaluationTapeInterface
      /// @{

      /// \copydoc codi::PositionalEvaluationTapeInterface::evaluate()
      /// <br> Implementation: Automatic adjoints management only involves bounds checking and resizing. Primal value
      /// tapes do not implement adjoints locking.
      CODI_INLINE void evaluate(Position const& start, Position const& end,
                                AdjointsManagement adjointsManagement = AdjointsManagement::Automatic) {
        if (AdjointsManagement::Automatic == adjointsManagement) {
          checkAdjointSize(indexManager.get().getLargestCreatedIndex());
        }

        codiAssert(indexManager.get().getLargestCreatedIndex() < (Identifier)adjoints.size());

        evaluate(start, end, adjoints.data());
      }

      /// \copydoc codi::PositionalEvaluationTapeInterface::resetTo(T_Position const&, bool, AdjointsManagement)
      /// <br> Implementation: Automatic adjoints management has no effect. Primal value tapes do not implement adjoints
      /// locking.
      CODI_INLINE void resetTo(Position const& pos, bool resetAdjoints = true,
                               AdjointsManagement adjointsManagement = AdjointsManagement::Automatic) {
        cast().internalResetPrimalValues(pos);

        Base::resetTo(pos, resetAdjoints, adjointsManagement);
      }

      /// @}
      /*******************************************************************************/
      /// @name Functions from PreaccumulationEvaluationTapeInterface
      /// @{

      /// \copydoc codi::PreaccumulationEvaluationTapeInterface::evaluateKeepState()
      /// <br> Implementation: Automatic adjoints management only involves bounds checking and resizing. Primal value
      /// tapes do not implement adjoints locking.
      void evaluateKeepState(Position const& start, Position const& end,
                             AdjointsManagement adjointsManagement = AdjointsManagement::Automatic) {
        if (AdjointsManagement::Automatic == adjointsManagement) {
          checkAdjointSize(indexManager.get().getLargestCreatedIndex());
        }

        codiAssert(indexManager.get().getLargestCreatedIndex() < (Identifier)adjoints.size());

        evaluateKeepState(start, end, adjoints.data());
      }

      /// \copydoc codi::PreaccumulationEvaluationTapeInterface::evaluateKeepState()
      template<typename AdjointVector>
      void evaluateKeepState(Position const& start, Position const& end, AdjointVector&& data) {
        internalEvaluateReverse<false>(start, end, std::forward<AdjointVector>(data));

        if (!TapeTypes::IsLinearIndexHandler) {
          evaluatePrimal(end, start);
        }
      }

      /// \copydoc codi::PreaccumulationEvaluationTapeInterface::evaluateForwardKeepState()
      /// <br> Implementation: Automatic adjoints management only involves bounds checking and resizing. Primal value
      /// tapes do not implement adjoints locking.
      void evaluateForwardKeepState(Position const& start, Position const& end,
                                    AdjointsManagement adjointsManagement = AdjointsManagement::Automatic) {
        if (AdjointsManagement::Automatic == adjointsManagement) {
          checkAdjointSize(indexManager.get().getLargestCreatedIndex());
        }

        codiAssert(indexManager.get().getLargestCreatedIndex() < (Identifier)adjoints.size());

        evaluateForwardKeepState(start, end, adjoints.data());
      }

      /// \copydoc codi::PreaccumulationEvaluationTapeInterface::evaluateForwardKeepState()
      template<typename AdjointVector>
      void evaluateForwardKeepState(Position const& start, Position const& end, AdjointVector&& data) {
        if (!TapeTypes::IsLinearIndexHandler) {
          cast().internalResetPrimalValues(end);
        }

        internalEvaluateForward<false>(start, end, std::forward<AdjointVector>(data));
      }

    protected:

      /******************************************************************************
       * Protected helper function for PrimalEvaluationTapeInterface
       */

      /// Wrapper helper for improved compiler optimizations.
      CODI_WRAP_FUNCTION(Wrap_internalEvaluatePrimal_EvalStatements, Impl::internalEvaluatePrimal_EvalStatements);

    public:

      /// @}
      /*******************************************************************************/
      /// @name Functions from PrimalEvaluationTapeInterface
      /// @{

      using Base::evaluatePrimal;

      /// \copydoc codi::PrimalEvaluationTapeInterface::evaluatePrimal()
      CODI_NO_INLINE void evaluatePrimal(Position const& start, Position const& end) {
        // TODO: implement primal value only accessor
        PrimalAdjointVectorAccess<Real, Identifier, Gradient*> primalAdjointAccess(adjoints.data(), primals.data());

        EventSystem<Impl>::notifyTapeEvaluateListeners(cast(), start, end, &primalAdjointAccess,
                                                       EventHints::EvaluationKind::Primal, EventHints::Endpoint::Begin);

        Wrap_internalEvaluatePrimal_EvalStatements evalFunc{};
        Base::llfByteData.evaluateForward(start, end, evalFunc, cast(), primals.data());

        EventSystem<Impl>::notifyTapeEvaluateListeners(cast(), start, end, &primalAdjointAccess,
                                                       EventHints::EvaluationKind::Primal, EventHints::Endpoint::End);
      }

      /// \copydoc codi::PrimalEvaluationTapeInterface::primal(T_Identifier const&)
      Real& primal(Identifier const& identifier) {
        return primals[identifier];
      }

      /// \copydoc codi::PrimalEvaluationTapeInterface::primal(T_Identifier const&) const
      Real const& primal(Identifier const& identifier) const {
        return primals[identifier];
      }

      /// \copydoc codi::PrimalEvaluationTapeInterface::getPrimalVector
      Real* getPrimalVector() {
        return primals.data();
      }

      /// @}
      /*******************************************************************************/
      /// @name Function from StatementEvaluatorInnerTapeInterface and StatementEvaluatorInnerTapeInterface
      /// @{

      /// Base class for statement call generators. Prepares the static construction context.
      template<typename Stmt, typename GenImpl>
      struct StatementCallGeneratorBase {
        public:
          using Lhs = typename Stmt::Lhs;                                     ///< \copydoc codi::AssignStatement::Lhs
          using Rhs = typename Stmt::Rhs;                                     ///< \copydoc codi::AssignStatement::Rhs
          using LhsReal = typename Lhs::Real;                                 ///< Primal value of lhs.
          using AggregateTraits = RealTraits::AggregatedTypeTraits<LhsReal>;  ///< Traits for aggregated type.

          using Constructor = ConstructStaticContextLogic<Rhs, Impl, 0, 0>;  ///< Static construction context.
          using StaticRhs = typename Constructor::ResultType;                ///< Static right hand side.

          template<size_t pos>
          using ExtractExpr =
              ArrayAccessExpression<LhsReal, pos, StaticRhs>;  ///< Extract expressions for aggregated types.

          /// Construct the statement in the static context.
          CODI_INLINE static StaticRhs constructStaticRhs(Real* CODI_RESTRICT primalVector,
                                                          PassiveReal const* CODI_RESTRICT const constantValues,
                                                          Identifier const* CODI_RESTRICT const identifiers) {
            return Constructor::construct(primalVector, identifiers, constantValues);
          }

          /// \copydoc codi::StatementEvaluatorTapeInterface::StatementCallGenerator::evaluate()
          template<typename... Args>
          CODI_INLINE static void internalEvaluate(Args&&... args) {
            size_t constexpr MaxActiveArgs = ExpressionTraits::NumberOfActiveTypeArguments<Rhs>::value;
            size_t constexpr MaxConstantArgs = ExpressionTraits::NumberOfConstantTypeArguments<Rhs>::value;
            size_t constexpr MaxOutputArgs = ExpressionTraits::NumberOfActiveTypeArguments<Lhs>::value;

            GenImpl::evaluateFull(GenImpl::evaluateInner, MaxOutputArgs, MaxActiveArgs, MaxConstantArgs,
                                  std::forward<Args>(args)...);
          }
      };

      /// \copydoc StatementEvaluatorTapeInterface::StatementCallGenerator
      template<StatementCall type, typename Stmt>
      struct StatementCallGenerator;

// Define macros for the definition of the evaluate functions in the Statement call generators.
#define STMT_COMMON_ARGS Config::ArgumentSize numberOfPassiveArguments, char *CODI_RESTRICT byteData
#define STMT_COMMON_CALL numberOfPassiveArguments, byteData

      /*******************************************************************************/
      /// WriteInformation implementation
      template<typename Stmt>
      struct StatementCallGenerator<StatementCall::WriteInformation, Stmt>
          : public StatementCallGeneratorBase<Stmt, StatementCallGenerator<StatementCall::WriteInformation, Stmt>> {
        public:
          /// Base class abbreviation.
          using Base = StatementCallGeneratorBase<Stmt, StatementCallGenerator<StatementCall::WriteInformation, Stmt>>;
          using StaticRhs = typename Base::StaticRhs;  ///< See StatementCallGeneratorBase.

          /// \copydoc codi::StatementEvaluatorInnerTapeInterface::StatementCallGenerator::evaluateInner()
          CODI_INLINE static void evaluateInner(size_t const& CODI_RESTRICT maxOutputArgs,
                                                size_t const& CODI_RESTRICT maxActiveArgs,
                                                size_t const& CODI_RESTRICT maxConstantArgs,
                                                WriteInfo& CODI_RESTRICT writeInfo, Real* CODI_RESTRICT primalVector,
                                                STMT_COMMON_ARGS) {
            StatementDataPointers pointers = {};
            pointers.populate(maxOutputArgs, maxActiveArgs, numberOfPassiveArguments, maxConstantArgs, byteData);

            if (Config::StatementInputTag != numberOfPassiveArguments) {
              for (Config::ArgumentSize curPos = 0; curPos < numberOfPassiveArguments; curPos += 1) {
                primalVector[curPos] = pointers.passiveValues[curPos];
              }
            }

            StaticRhs staticRhs =
                Base::constructStaticRhs(primalVector, pointers.constantValues, pointers.rhsIdentifiers);

            writeInfo.numberOfOutputArguments = maxOutputArgs;
            writeInfo.numberOfActiveArguments = maxActiveArgs;
            writeInfo.numberOfConstantArguments = maxConstantArgs;

            // The Impl template name is added to simplify the writer and reader interfaces EvalHandles. It represents
            // the first and second template args in the createHandle method. For the manual push, the second Impl
            // template is instead replaced with JacobianGenerator<size>.

            writeInfo.stmtExpression = "Impl, Impl, ";
            writeInfo.stmtExpression += demangleName<Stmt>();
            MathStatementGenLogic<Identifier> mathGen(Config::MaxArgumentSize);
            mathGen.eval(staticRhs, writeInfo.mathRepresentation);
          }

          /// \copydoc codi::StatementEvaluatorInnerTapeInterface::StatementCallGenerator::evaluateFull()
          template<typename Func, typename... Args>
          CODI_INLINE static void evaluateFull(Func const& func, Args&&... args) {
            func(std::forward<Args>(args)...);
          }

          /// \copydoc codi::StatementEvaluatorTapeInterface::StatementCallGenerator::evaluate()
          CODI_INLINE static void evaluate(WriteInfo& CODI_RESTRICT writeInfo, Real* CODI_RESTRICT primalVector,
                                           STMT_COMMON_ARGS) {
            Base::internalEvaluate(writeInfo, primalVector, STMT_COMMON_CALL);
          }
      };

      /*******************************************************************************/
      /// IterateInputs implementation
      template<typename Stmt>
      struct StatementCallGenerator<StatementCall::IterateInputs, Stmt>
          : public StatementCallGeneratorBase<Stmt, StatementCallGenerator<StatementCall::IterateInputs, Stmt>> {
        public:
          /// Base class abbreviation.
          using Base = StatementCallGeneratorBase<Stmt, StatementCallGenerator<StatementCall::IterateInputs, Stmt>>;
          using StaticRhs = typename Base::StaticRhs;  ///< See StatementCallGeneratorBase.

          /// See LowLevelFunctionEntry::IterCallback.
          using IterCallback = typename LowLevelFunctionEntry<Impl, Real, Identifier>::IterCallback;

          /// \copydoc codi::StatementEvaluatorInnerTapeInterface::StatementCallGenerator::evaluateInner()
          CODI_INLINE static void evaluateInner(size_t const& CODI_RESTRICT maxOutputArgs,
                                                size_t const& CODI_RESTRICT maxActiveArgs,
                                                size_t const& CODI_RESTRICT maxConstantArgs,
                                                size_t& CODI_RESTRICT linearAdjointPos, IterCallback func,
                                                void* userData, STMT_COMMON_ARGS) {
            CODI_UNUSED(linearAdjointPos);

            StatementDataPointers pointers = {};
            pointers.populate(maxOutputArgs, maxActiveArgs, numberOfPassiveArguments, maxConstantArgs, byteData);

            for (Config::ArgumentSize i = 0; i < maxActiveArgs; i += 1) {
              if (pointers.rhsIdentifiers[i] >= (Identifier)Config::MaxArgumentSize) {
                func(&pointers.rhsIdentifiers[i], userData);
              }
            }
          }

          /// \copydoc codi::StatementEvaluatorInnerTapeInterface::StatementCallGenerator::evaluateFull()
          template<typename Func, typename... Args>
          CODI_INLINE static void evaluateFull(Func const& func, Args&&... args) {
            func(std::forward<Args>(args)...);
          }

          /// \copydoc codi::StatementEvaluatorTapeInterface::StatementCallGenerator::evaluate()
          CODI_INLINE static void evaluate(size_t& CODI_RESTRICT linearAdjointPos, IterCallback func, void* userData,
                                           STMT_COMMON_ARGS) {
            Base::internalEvaluate(linearAdjointPos, func, userData, STMT_COMMON_CALL);
          }
      };

      /*******************************************************************************/
      /// IterateOutputs implementation
      template<typename Stmt>
      struct StatementCallGenerator<StatementCall::IterateOutputs, Stmt>
          : public StatementCallGeneratorBase<Stmt, StatementCallGenerator<StatementCall::IterateOutputs, Stmt>> {
        public:
          /// Base class abbreviation.
          using Base = StatementCallGeneratorBase<Stmt, StatementCallGenerator<StatementCall::IterateOutputs, Stmt>>;
          using StaticRhs = typename Base::StaticRhs;  ///< See StatementCallGeneratorBase.

          /// See LowLevelFunctionEntry::IterCallback.
          using IterCallback = typename LowLevelFunctionEntry<Impl, Real, Identifier>::IterCallback;

          /// \copydoc codi::StatementEvaluatorInnerTapeInterface::StatementCallGenerator::evaluateInner()
          CODI_INLINE static void evaluateInner(size_t const& CODI_RESTRICT maxOutputArgs,
                                                size_t const& CODI_RESTRICT maxActiveArgs,
                                                size_t const& CODI_RESTRICT maxConstantArgs,
                                                size_t& CODI_RESTRICT linearAdjointPos, IterCallback func,
                                                void* userData, STMT_COMMON_ARGS) {
            StatementDataPointers pointers = {};
            pointers.populate(maxOutputArgs, maxActiveArgs, numberOfPassiveArguments, maxConstantArgs, byteData);

            for (size_t iLhs = 0; iLhs < maxOutputArgs; iLhs += 1) {
              if (LinearIndexHandling) {
                Identifier lhsIdentifier = linearAdjointPos + 1 + iLhs;
                func(&lhsIdentifier, userData);

                codiAssert(lhsIdentifier == (Identifier)(linearAdjointPos + 1 + iLhs));
              } else {
                func(&pointers.lhsIdentifiers[iLhs], userData);
              }
            }
          }

          /// \copydoc codi::StatementEvaluatorInnerTapeInterface::StatementCallGenerator::evaluateFull()
          template<typename Func, typename... Args>
          CODI_INLINE static void evaluateFull(Func const& func, Args&&... args) {
            func(std::forward<Args>(args)...);
          }

          /// \copydoc codi::StatementEvaluatorTapeInterface::StatementCallGenerator::evaluate()
          CODI_INLINE static void evaluate(size_t& CODI_RESTRICT linearAdjointPos, IterCallback func, void* userData,
                                           STMT_COMMON_ARGS) {
            Base::internalEvaluate(linearAdjointPos, func, userData, STMT_COMMON_CALL);
          }
      };

      /*******************************************************************************/
      /// ClearAdjoints implementation
      template<typename Stmt>
      struct StatementCallGenerator<StatementCall::ClearAdjoints, Stmt>
          : public StatementCallGeneratorBase<Stmt, StatementCallGenerator<StatementCall::ClearAdjoints, Stmt>> {
        public:
          /// Base class abbreviation.
          using Base = StatementCallGeneratorBase<Stmt, StatementCallGenerator<StatementCall::ClearAdjoints, Stmt>>;
          using StaticRhs = typename Base::StaticRhs;  ///< See StatementCallGeneratorBase.

          /// \copydoc codi::StatementEvaluatorInnerTapeInterface::StatementCallGenerator::evaluateInner()
          CODI_INLINE static void evaluateInner() {}

          /// \copydoc codi::StatementEvaluatorInnerTapeInterface::StatementCallGenerator::evaluateFull()
          template<typename Func>
          CODI_INLINE static void evaluateFull(Func const& evalInner, size_t const& CODI_RESTRICT maxOutputArgs,
                                               size_t const& CODI_RESTRICT maxActiveArgs,
                                               size_t const& CODI_RESTRICT maxConstantArgs,
                                               ADJOINT_VECTOR_TYPE* CODI_RESTRICT adjointVector, STMT_COMMON_ARGS) {
            CODI_UNUSED(evalInner);

            StatementDataPointers pointers = {};
            pointers.populate(maxOutputArgs, maxActiveArgs, numberOfPassiveArguments, maxConstantArgs, byteData);

            for (size_t iLhs = 0; iLhs < maxOutputArgs; iLhs += 1) {
              // Only called from reuse tape.
              Identifier lhsIdentifier = pointers.lhsIdentifiers[iLhs];

#if CODI_VariableAdjointInterfaceInPrimalTapes
              adjointVector->resetAdjointVec(lhsIdentifier);
#else
              adjointVector[lhsIdentifier] = Gradient();
#endif
            }
          }

          /// \copydoc codi::StatementEvaluatorTapeInterface::StatementCallGenerator::evaluate()
          CODI_INLINE static void evaluate(ADJOINT_VECTOR_TYPE* CODI_RESTRICT adjointVector, STMT_COMMON_ARGS) {
            Base::internalEvaluate(adjointVector, STMT_COMMON_CALL);
          }
      };

      /*******************************************************************************/
      /// Forward implementation
      template<typename Stmt>
      struct StatementCallGenerator<StatementCall::Forward, Stmt>
          : public StatementCallGeneratorBase<Stmt, StatementCallGenerator<StatementCall::Forward, Stmt>> {
        public:
          /// Base class abbreviation.
          using Base = StatementCallGeneratorBase<Stmt, StatementCallGenerator<StatementCall::Forward, Stmt>>;
          using StaticRhs = typename Base::StaticRhs;  ///< See StatementCallGeneratorBase.

          template<size_t pos>
          using ExtractExpr = typename Base::template ExtractExpr<pos>;  ///< See StatementCallGeneratorBase.

          /// \copydoc codi::StatementEvaluatorInnerTapeInterface::StatementCallGenerator::evaluateInner()
          CODI_INLINE static void evaluateInner(Real* CODI_RESTRICT primalVector,
                                                ADJOINT_VECTOR_TYPE* CODI_RESTRICT adjointVector,
                                                Real* CODI_RESTRICT lhsPrimals, Gradient* lhsTangents,
                                                PassiveReal const* CODI_RESTRICT const constantValues,
                                                Identifier const* CODI_RESTRICT const rhsIdentifiers) {
            StaticRhs staticRhs = Base::constructStaticRhs(primalVector, constantValues, rhsIdentifiers);

            IncrementForwardLogic incrementForward;

            static_for<Base::AggregateTraits::Elements>([&](auto i) CODI_LAMBDA_INLINE {
#if CODI_VariableAdjointInterfaceInPrimalTapes
              adjointVector->setActiveVariableForIndirectAccess(i.value);
#else
              lhsTangents[i.value] = Gradient();
#endif
              ExtractExpr<i.value> expr(staticRhs);

              incrementForward.eval(expr, Real(1.0), lhsTangents[i.value], adjointVector);
              lhsPrimals[i.value] = expr.getValue();
            });
          }

          /// \copydoc codi::StatementEvaluatorInnerTapeInterface::StatementCallGenerator::evaluateFull()
          template<typename Func>
          CODI_INLINE static void evaluateFull(Func const& evalInner, size_t const& CODI_RESTRICT maxOutputArgs,
                                               size_t const& CODI_RESTRICT maxActiveArgs,
                                               size_t const& CODI_RESTRICT maxConstantArgs, Impl& CODI_RESTRICT tape,
                                               Real* CODI_RESTRICT lhsPrimals, Gradient* CODI_RESTRICT lhsTangents,
                                               Real* CODI_RESTRICT primalVector,
                                               ADJOINT_VECTOR_TYPE* CODI_RESTRICT adjointVector,
                                               size_t& CODI_RESTRICT linearAdjointPos, STMT_COMMON_ARGS) {
            StatementDataPointers pointers = {};
            pointers.populate(maxOutputArgs, maxActiveArgs, numberOfPassiveArguments, maxConstantArgs, byteData);

            for (Config::ArgumentSize curPos = 0; curPos < numberOfPassiveArguments; curPos += 1) {
              primalVector[curPos] = pointers.passiveValues[curPos];
            }

            evalInner(primalVector, adjointVector, lhsPrimals, lhsTangents, pointers.constantValues,
                      pointers.rhsIdentifiers);

            for (size_t iLhs = 0; iLhs < maxOutputArgs; iLhs += 1) {
              Identifier lhsIdentifier;
              if (LinearIndexHandling) {
                linearAdjointPos += 1;
                lhsIdentifier = linearAdjointPos;
              } else {
                lhsIdentifier = pointers.lhsIdentifiers[iLhs];
              }

              if (!LinearIndexHandling) {
                pointers.oldLhsValues[iLhs] = primalVector[lhsIdentifier];
              }

              primalVector[lhsIdentifier] = lhsPrimals[iLhs];

#if CODI_VariableAdjointInterfaceInPrimalTapes
              adjointVector->setLhsTangent(lhsIdentifier);
              EventSystem<Impl>::notifyStatementEvaluateListeners(tape, lhsIdentifier, adjointVector->getVectorSize(),
                                                                  adjointVector->getAdjointVec(lhsIdentifier));
#else
              adjointVector[lhsIdentifier] = lhsTangents[iLhs];
              EventSystem<Impl>::notifyStatementEvaluateListeners(tape, lhsIdentifier, GradientTraits::dim<Gradient>(),
                                                                  GradientTraits::toArray(lhsTangents[iLhs]).data());
#endif
              EventSystem<Impl>::notifyStatementEvaluatePrimalListeners(tape, lhsIdentifier,
                                                                        primalVector[lhsIdentifier]);
            }
          }

          /// \copydoc codi::StatementEvaluatorTapeInterface::StatementCallGenerator::evaluate()
          CODI_INLINE static void evaluate(Impl& CODI_RESTRICT tape, Real* CODI_RESTRICT lhsPrimals,
                                           Gradient* CODI_RESTRICT lhsTangents, Real* CODI_RESTRICT primalVector,
                                           ADJOINT_VECTOR_TYPE* CODI_RESTRICT adjointVector,
                                           size_t& CODI_RESTRICT linearAdjointPos, STMT_COMMON_ARGS) {
            Base::internalEvaluate(tape, lhsPrimals, lhsTangents, primalVector, adjointVector, linearAdjointPos,
                                   STMT_COMMON_CALL);
          }
      };

      /*******************************************************************************/
      /// Primal implementation
      template<typename Stmt>
      struct StatementCallGenerator<StatementCall::Primal, Stmt>
          : public StatementCallGeneratorBase<Stmt, StatementCallGenerator<StatementCall::Primal, Stmt>> {
        public:
          /// Base class abbreviation.
          using Base = StatementCallGeneratorBase<Stmt, StatementCallGenerator<StatementCall::Primal, Stmt>>;
          using StaticRhs = typename Base::StaticRhs;  ///< See StatementCallGeneratorBase.

          template<size_t pos>
          using ExtractExpr = typename Base::template ExtractExpr<pos>;  ///< See StatementCallGeneratorBase.

          /// \copydoc codi::StatementEvaluatorInnerTapeInterface::StatementCallGenerator::evaluateInner()
          CODI_INLINE static void evaluateInner(Real* CODI_RESTRICT lhsPrimals, Real* CODI_RESTRICT primalVector,
                                                PassiveReal const* CODI_RESTRICT const constantValues,
                                                Identifier const* CODI_RESTRICT const rhsIdentifiers) {
            StaticRhs staticRhs = Base::constructStaticRhs(primalVector, constantValues, rhsIdentifiers);

            static_for<Base::AggregateTraits::Elements>([&](auto i) CODI_LAMBDA_INLINE {
              ExtractExpr<i.value> expr(staticRhs);

              lhsPrimals[i.value] = expr.getValue();
            });
          }

          /// \copydoc codi::StatementEvaluatorInnerTapeInterface::StatementCallGenerator::evaluateFull()
          template<typename Func>
          CODI_INLINE static void evaluateFull(Func const& evalInner, size_t const& CODI_RESTRICT maxOutputArgs,
                                               size_t const& CODI_RESTRICT maxActiveArgs,
                                               size_t const& CODI_RESTRICT maxConstantArgs, Impl& CODI_RESTRICT tape,
                                               Real* CODI_RESTRICT lhsPrimals, Real* CODI_RESTRICT primalVector,
                                               size_t& CODI_RESTRICT linearAdjointPos, STMT_COMMON_ARGS) {
            StatementDataPointers pointers = {};
            pointers.populate(maxOutputArgs, maxActiveArgs, numberOfPassiveArguments, maxConstantArgs, byteData);

            for (Config::ArgumentSize curPos = 0; curPos < numberOfPassiveArguments; curPos += 1) {
              primalVector[curPos] = pointers.passiveValues[curPos];
            }

            evalInner(lhsPrimals, primalVector, pointers.constantValues, pointers.rhsIdentifiers);

            for (size_t iLhs = 0; iLhs < maxOutputArgs; iLhs += 1) {
              Identifier lhsIdentifier;
              if (LinearIndexHandling) {
                linearAdjointPos += 1;
                lhsIdentifier = linearAdjointPos;
              } else {
                lhsIdentifier = pointers.lhsIdentifiers[iLhs];
              }

              if (!LinearIndexHandling) {
                pointers.oldLhsValues[iLhs] = primalVector[lhsIdentifier];
              }

              primalVector[lhsIdentifier] = lhsPrimals[iLhs];

              EventSystem<Impl>::notifyStatementEvaluatePrimalListeners(tape, lhsIdentifier,
                                                                        primalVector[lhsIdentifier]);
            }
          }

          /// \copydoc codi::StatementEvaluatorTapeInterface::StatementCallGenerator::evaluate()
          CODI_INLINE static void evaluate(Impl& CODI_RESTRICT tape, Real* CODI_RESTRICT lhsPrimals,
                                           Real* CODI_RESTRICT primalVector, size_t& CODI_RESTRICT linearAdjointPos,
                                           STMT_COMMON_ARGS) {
            Base::internalEvaluate(tape, lhsPrimals, primalVector, linearAdjointPos, STMT_COMMON_CALL);
          }
      };

      /*******************************************************************************/
      /// Reverse implementation
      template<typename Stmt>
      struct StatementCallGenerator<StatementCall::Reverse, Stmt>
          : public StatementCallGeneratorBase<Stmt, StatementCallGenerator<StatementCall::Reverse, Stmt>> {
        public:
          /// Base class abbreviation.
          using Base = StatementCallGeneratorBase<Stmt, StatementCallGenerator<StatementCall::Reverse, Stmt>>;
          using StaticRhs = typename Base::StaticRhs;  ///< See StatementCallGeneratorBase.

          template<size_t pos>
          using ExtractExpr = typename Base::template ExtractExpr<pos>;  ///< See StatementCallGeneratorBase.

          /// \copydoc codi::StatementEvaluatorInnerTapeInterface::StatementCallGenerator::evaluateInner()
          CODI_INLINE static void evaluateInner(Real* CODI_RESTRICT primalVector,
                                                ADJOINT_VECTOR_TYPE* CODI_RESTRICT adjointVector,
                                                Gradient* CODI_RESTRICT lhsAdjoints,
                                                PassiveReal const* CODI_RESTRICT const constantValues,
                                                Identifier const* CODI_RESTRICT const rhsIdentifiers) {
            StaticRhs staticRhs = Base::constructStaticRhs(primalVector, constantValues, rhsIdentifiers);

            IncrementReversalLogic incrementReverse;
            static_for<Base::AggregateTraits::Elements>([&](auto i) CODI_LAMBDA_INLINE {
#if CODI_VariableAdjointInterfaceInPrimalTapes
              adjointVector->setActiveVariableForIndirectAccess(i.value);
#endif
              ExtractExpr<i.value> expr(staticRhs);
              incrementReverse.eval(expr, Real(1.0), const_cast<Gradient const&>(lhsAdjoints[i.value]), adjointVector);
            });
          }

          /// \copydoc codi::StatementEvaluatorInnerTapeInterface::StatementCallGenerator::evaluateFull()
          template<typename Func>
          CODI_INLINE static void evaluateFull(Func const& evalInner, size_t const& CODI_RESTRICT maxOutputArgs,
                                               size_t const& CODI_RESTRICT maxActiveArgs,
                                               size_t const& CODI_RESTRICT maxConstantArgs, Impl& CODI_RESTRICT tape,
                                               Gradient* CODI_RESTRICT lhsAdjoints, Real* CODI_RESTRICT primalVector,
                                               ADJOINT_VECTOR_TYPE* CODI_RESTRICT adjointVector,
                                               size_t& CODI_RESTRICT linearAdjointPos, STMT_COMMON_ARGS) {
            StatementDataPointers pointers = {};
            pointers.populate(maxOutputArgs, maxActiveArgs, numberOfPassiveArguments, maxConstantArgs, byteData);

            if (LinearIndexHandling) {
              linearAdjointPos -= maxOutputArgs;
            }

            bool allZero = true;
            for (size_t iLhs = 0; iLhs < maxOutputArgs; iLhs += 1) {
              Identifier lhsIdentifier;
              if (LinearIndexHandling) {
                lhsIdentifier = linearAdjointPos + 1 + iLhs;
              } else {
                lhsIdentifier = pointers.lhsIdentifiers[iLhs];
              }

#if CODI_VariableAdjointInterfaceInPrimalTapes
              EventSystem<Impl>::notifyStatementEvaluateListeners(tape, lhsIdentifier, adjointVector->getVectorSize(),
                                                                  adjointVector->getAdjointVec(lhsIdentifier));
              adjointVector->setActiveVariableForIndirectAccess(iLhs);
              adjointVector->setLhsAdjoint(lhsIdentifier);
              allZero &= adjointVector->isLhsZero();
#else
              lhsAdjoints[iLhs] = adjointVector[lhsIdentifier];
              EventSystem<Impl>::notifyStatementEvaluateListeners(tape, lhsIdentifier, GradientTraits::dim<Gradient>(),
                                                                  GradientTraits::toArray(lhsAdjoints[iLhs]).data());
              adjointVector[lhsIdentifier] = Gradient();
              allZero &= RealTraits::isTotalZero(lhsAdjoints[iLhs]);
#endif

              EventSystem<Impl>::notifyStatementEvaluatePrimalListeners(tape, lhsIdentifier,
                                                                        primalVector[lhsIdentifier]);
              if (!LinearIndexHandling) {
                primalVector[lhsIdentifier] = pointers.oldLhsValues[iLhs];
              }
            }

            if (CODI_ENABLE_CHECK(Config::SkipZeroAdjointEvaluation, !allZero)) {
              for (Config::ArgumentSize curPos = 0; curPos < numberOfPassiveArguments; curPos += 1) {
                primalVector[curPos] = pointers.passiveValues[curPos];
              }

              evalInner(primalVector, adjointVector, lhsAdjoints, pointers.constantValues, pointers.rhsIdentifiers);
            }
          }

          /// \copydoc codi::StatementEvaluatorTapeInterface::StatementCallGenerator::evaluate()
          CODI_INLINE static void evaluate(Impl& CODI_RESTRICT tape, Gradient* CODI_RESTRICT lhsAdjoints,
                                           Real* CODI_RESTRICT primalVector,
                                           ADJOINT_VECTOR_TYPE* CODI_RESTRICT adjointVector,
                                           size_t& CODI_RESTRICT linearAdjointPos, STMT_COMMON_ARGS) {
            Base::internalEvaluate(tape, lhsAdjoints, primalVector, adjointVector, linearAdjointPos, STMT_COMMON_CALL);
          }
      };

      /*******************************************************************************/
      /// ResetPrimal implementation
      template<typename Stmt>
      struct StatementCallGenerator<StatementCall::ResetPrimals, Stmt>
          : public StatementCallGeneratorBase<Stmt, StatementCallGenerator<StatementCall::ResetPrimals, Stmt>> {
        public:
          /// Base class abbreviation.
          using Base = StatementCallGeneratorBase<Stmt, StatementCallGenerator<StatementCall::ResetPrimals, Stmt>>;
          using StaticRhs = typename Base::StaticRhs;  ///< See StatementCallGeneratorBase.

          /// \copydoc codi::StatementEvaluatorInnerTapeInterface::StatementCallGenerator::evaluateInner()
          CODI_INLINE static void evaluateInner() {}

          /// \copydoc codi::StatementEvaluatorInnerTapeInterface::StatementCallGenerator::evaluateFull()
          template<typename Func>
          CODI_INLINE static void evaluateFull(Func const& evalInner, size_t const& CODI_RESTRICT maxOutputArgs,
                                               size_t const& CODI_RESTRICT maxActiveArgs,
                                               size_t const& CODI_RESTRICT maxConstantArgs,
                                               Real* CODI_RESTRICT primalVector, STMT_COMMON_ARGS) {
            CODI_UNUSED(evalInner);

            StatementDataPointers pointers = {};
            pointers.populate(maxOutputArgs, maxActiveArgs, numberOfPassiveArguments, maxConstantArgs, byteData);

            // Only called from reuse tape.
            for (size_t iLhs = 0; iLhs < maxOutputArgs; iLhs += 1) {
              Identifier lhsIdentifier = pointers.lhsIdentifiers[iLhs];

              primalVector[lhsIdentifier] = pointers.oldLhsValues[iLhs];
            }
          }

          /// \copydoc codi::StatementEvaluatorTapeInterface::StatementCallGenerator::evaluate()
          CODI_INLINE static void evaluate(Real* CODI_RESTRICT primalVector, STMT_COMMON_ARGS) {
            Base::internalEvaluate(primalVector, STMT_COMMON_CALL);
          }
      };

      /// @}

    private:

      CODI_INLINE void checkAdjointSize(Identifier const& identifier) {
        if (identifier >= (Identifier)adjoints.size()) {
          resizeAdjointsVector();
        }
      }

      CODI_INLINE void checkPrimalSize(bool generatedNewIndex) {
        if (generatedNewIndex && indexManager.get().getLargestCreatedIndex() >= (Identifier)primals.size()) {
          resizePrimalVector();
        }
      }

      CODI_NO_INLINE void resizeAdjointsVector() {
        // overallocate as next multiple of Config::ChunkSize
        adjoints.resize(getNextMultiple((size_t)indexManager.get().getLargestCreatedIndex() + 1, Config::ChunkSize));
      }

      CODI_NO_INLINE void resizePrimalVector() {
        // overallocate as next multiple of Config::ChunkSize
        primals.resize(getNextMultiple((size_t)indexManager.get().getLargestCreatedIndex() + 1, Config::ChunkSize));
      }
  };

  /// Specialized for NumberOfActiveTypeArguments and NumberOfConstantTypeArguments.
  template<size_t size>
  struct JacobianExpression {};

  /// Specialization for manual statement pushes of the used expression type.
  template<size_t size>
  struct ExpressionTraits::NumberOfActiveTypeArguments<JacobianExpression<size>> {
      static size_t constexpr value = size;  ///< Number of arguments.
  };

  /// Specialization for manual statement pushes of the used expression type.
  template<size_t size>
  struct ExpressionTraits::NumberOfConstantTypeArguments<JacobianExpression<size>> {
      static size_t constexpr value = 0;  ///< Always zero.
  };

  /// Implements StatementEvaluatorTapeInterface and StatementEvaluatorInnerTapeInterface
  /// @tparam T_size  Number of arguments.
  template<typename T_TapeImpl, size_t T_size>
  struct JacobianStatementGenerator : public StatementEvaluatorTapeInterface,
                                      public StatementEvaluatorInnerTapeInterface {
    public:

      using TapeImpl = CODI_DD(
          T_TapeImpl, CODI_T(PrimalValueBaseTape<PrimalValueTapeTypes<double, double, IndexManagerInterface<int>,
                                                                      StatementEvaluatorInterface, DefaultChunkedData>,
                                                 void>));  ///< PrimalValueBaseTape.
      static size_t constexpr size = T_size;               ///< See JacobianStatementGenerator

      using Real = typename TapeImpl::Real;                ///< See PrimalValueBaseTape.
      using Gradient = typename TapeImpl::Gradient;        ///< See PrimalValueBaseTape.
      using Identifier = typename TapeImpl::Identifier;    ///< See PrimalValueBaseTape.
      using PassiveReal = typename TapeImpl::PassiveReal;  ///< See PrimalValueBaseTape.

      using StatementDataPointers = typename TapeImpl::StatementDataPointers;  ///< Defined in PrimalValueBaseTape.

      static size_t constexpr LinearIndexHandling = TapeImpl::LinearIndexHandling;  ///< See PrimalValueBaseTape.

      /// Base for statement generators.
      template<typename Rhs, typename Impl>
      using StatementCallGeneratorBase = typename TapeImpl::template StatementCallGeneratorBase<Rhs, Impl>;

      /*******************************************************************************/
      /// @name Implementation of StatementEvaluatorTapeInterface and StatementEvaluatorInnerTapeInterface
      /// @{

      /// \copydoc StatementEvaluatorTapeInterface::StatementCallGenerator
      template<StatementCall type, typename Stmt>
      struct StatementCallGenerator;

      /// WriteInformation implementation
      template<typename Stmt>
      struct StatementCallGenerator<StatementCall::WriteInformation, Stmt>
          : public StatementCallGeneratorBase<Stmt, StatementCallGenerator<StatementCall::WriteInformation, Stmt>> {
        public:
          /// Base class abbreviation.
          using Base = StatementCallGeneratorBase<Stmt, StatementCallGenerator<StatementCall::WriteInformation, Stmt>>;

          /// \copydoc codi::StatementEvaluatorInnerTapeInterface::StatementCallGenerator::evaluateInner()
          CODI_INLINE static void evaluateInner() {}

          /// \copydoc codi::StatementEvaluatorInnerTapeInterface::StatementCallGenerator::evaluateFull()
          template<typename Func>
          CODI_INLINE static void evaluateFull(Func const& func, size_t const& maxOutputArgs,
                                               size_t const& maxActiveArgs, size_t const& maxConstantArgs,
                                               WriteInfo& writeInfo, Real* primalVector, STMT_COMMON_ARGS) {
            CODI_UNUSED(func, primalVector, numberOfPassiveArguments, byteData);

            writeInfo.numberOfOutputArguments = maxOutputArgs;
            writeInfo.numberOfActiveArguments = maxActiveArgs;
            writeInfo.numberOfConstantArguments = maxConstantArgs;
            writeInfo.stmtExpression = "Impl, typename Impl::template JacobianStatementGenerator<" +
                                       std::to_string(size) + ">, codi::JacobianExpression<" + std::to_string(size) +
                                       ">";
            writeInfo.mathRepresentation = "Jacobian statement";
          }

          /// \copydoc codi::StatementEvaluatorTapeInterface::StatementCallGenerator::evaluate()
          CODI_INLINE static void evaluate(WriteInfo& CODI_RESTRICT writeInfo, Real* CODI_RESTRICT primalVector,
                                           STMT_COMMON_ARGS) {
            Base::internalEvaluate(writeInfo, primalVector, STMT_COMMON_CALL);
          }
      };

      /// WriteInformation implementation
      template<typename Stmt>
      struct StatementCallGenerator<StatementCall::ClearAdjoints, Stmt>
          : public StatementCallGeneratorBase<Stmt, StatementCallGenerator<StatementCall::ClearAdjoints, Stmt>> {
        public:
          /// Base class abbreviation.
          using Base = StatementCallGeneratorBase<Stmt, StatementCallGenerator<StatementCall::ClearAdjoints, Stmt>>;

          /// \copydoc codi::StatementEvaluatorInnerTapeInterface::StatementCallGenerator::evaluateInner()
          CODI_INLINE static void evaluateInner() {}

          /// \copydoc codi::StatementEvaluatorInnerTapeInterface::StatementCallGenerator::evaluateFull()
          template<typename Func>
          CODI_INLINE static void evaluateFull(Func const& func, size_t const& maxOutputArgs,
                                               size_t const& maxActiveArgs, size_t const& maxConstantArgs,
                                               ADJOINT_VECTOR_TYPE* adjointVector, STMT_COMMON_ARGS) {
            CODI_UNUSED(func);

            StatementDataPointers pointers = {};
            pointers.populate(maxOutputArgs, maxActiveArgs, numberOfPassiveArguments, maxConstantArgs, byteData);

            codiAssert(1 == maxOutputArgs);

            // Only called from reuse tape.
            Identifier lhsIdentifier = pointers.lhsIdentifiers[0];

#if CODI_VariableAdjointInterfaceInPrimalTapes
            adjointVector->resetAdjointVec(lhsIdentifier);
#else
            adjointVector[lhsIdentifier] = Gradient();
#endif
          }

          /// \copydoc codi::StatementEvaluatorTapeInterface::StatementCallGenerator::evaluate()
          CODI_INLINE static void evaluate(ADJOINT_VECTOR_TYPE* CODI_RESTRICT adjointVector, STMT_COMMON_ARGS) {
            Base::internalEvaluate(adjointVector, STMT_COMMON_CALL);
          }
      };

      /// Forward implementation
      template<typename Stmt>
      struct StatementCallGenerator<StatementCall::Forward, Stmt> {
          /// Throws exception.
          static void evaluateInner() {
            CODI_EXCEPTION("Forward evaluation of jacobian statement not possible.");
          }

          /// Throws exception.
          static void evaluate() {
            CODI_EXCEPTION("Forward evaluation of jacobian statement not possible.");
          }
      };

      /// Primal implementation
      template<typename Stmt>
      struct StatementCallGenerator<StatementCall::Primal, Stmt> {
          /// Throws exception.
          static void evaluateInner() {
            CODI_EXCEPTION("Primal evaluation of jacobian statement not possible.");
          }

          /// Throws exception.
          static void evaluate() {
            CODI_EXCEPTION("Primal evaluation of jacobian statement not possible.");
          }
      };

      /// Reverse implementation
      template<typename Stmt>
      struct StatementCallGenerator<StatementCall::Reverse, Stmt>
          : public StatementCallGeneratorBase<Stmt, StatementCallGenerator<StatementCall::Reverse, Stmt>> {
        public:
          /// Base class abbreviation.
          using Base = StatementCallGeneratorBase<Stmt, StatementCallGenerator<StatementCall::Reverse, Stmt>>;

          /// \copydoc codi::StatementEvaluatorInnerTapeInterface::StatementCallGenerator::evaluateInner
          static void evaluateInner(Real* primalVector, ADJOINT_VECTOR_TYPE* adjointVector, Gradient* lhsAdjoints,
                                    PassiveReal const* const constantValues, Identifier const* const rhsIdentifiers) {
            CODI_UNUSED(constantValues);

            evalJacobianReverse(adjointVector, lhsAdjoints[0], primalVector, rhsIdentifiers);
          }

          /// \copydoc codi::StatementEvaluatorTapeInterface::StatementCallGenerator::evaluate
          static void evaluate(TapeImpl& tape, Gradient* lhsAdjoints, Real* primalVector,
                               ADJOINT_VECTOR_TYPE* adjointVector, size_t& linearAdjointPos, STMT_COMMON_ARGS) {
            CODI_UNUSED(primalVector, lhsAdjoints, numberOfPassiveArguments);

            StatementDataPointers pointers = {};
            pointers.populate(1, size, size, 0, byteData);

            Identifier lhsIdentifier;
            if (LinearIndexHandling) {
              lhsIdentifier = linearAdjointPos;
            } else {
              lhsIdentifier = pointers.lhsIdentifiers[0];
            }

#if CODI_VariableAdjointInterfaceInPrimalTapes
            Gradient const lhsAdjoint = {};
            EventSystem<TapeImpl>::notifyStatementEvaluateListeners(tape, lhsIdentifier, adjointVector->getVectorSize(),
                                                                    adjointVector->getAdjointVec(lhsIdentifier));
            adjointVector->setLhsAdjoint(lhsIdentifier);
#else
            Gradient const lhsAdjoint = adjointVector[lhsIdentifier];
            EventSystem<TapeImpl>::notifyStatementEvaluateListeners(
                tape, lhsIdentifier, GradientTraits::dim<Gradient>(), GradientTraits::toArray(lhsAdjoint).data());
            adjointVector[lhsIdentifier] = Gradient();
#endif

            EventSystem<TapeImpl>::notifyStatementEvaluatePrimalListeners(tape, lhsIdentifier,
                                                                          primalVector[lhsIdentifier]);
            if (!LinearIndexHandling) {
              primalVector[lhsIdentifier] = pointers.oldLhsValues[0];
            }

            evalJacobianReverse(adjointVector, lhsAdjoint, pointers.passiveValues, pointers.rhsIdentifiers);
          }
      };

      /// ResetPrimal implementation
      template<typename Stmt>
      struct StatementCallGenerator<StatementCall::ResetPrimals, Stmt>
          : public StatementCallGeneratorBase<Stmt, StatementCallGenerator<StatementCall::ResetPrimals, Stmt>> {
        public:
          /// Base class abbreviation.
          using Base = StatementCallGeneratorBase<Stmt, StatementCallGenerator<StatementCall::ResetPrimals, Stmt>>;

          /// \copydoc codi::StatementEvaluatorInnerTapeInterface::StatementCallGenerator::evaluateInner()
          CODI_INLINE static void evaluateInner() {}

          /// \copydoc codi::StatementEvaluatorInnerTapeInterface::StatementCallGenerator::evaluateFull()
          template<typename Func>
          CODI_INLINE static void evaluateFull(Func const& evalInner, size_t const& maxOutputArgs,
                                               size_t const& maxActiveArgs, size_t const& maxConstantArgs,
                                               Real* primalVector, STMT_COMMON_ARGS) {
            CODI_UNUSED(evalInner, maxOutputArgs, maxActiveArgs, maxConstantArgs, numberOfPassiveArguments);

            StatementDataPointers pointers = {};
            pointers.populate(1, size, size, 0, byteData);

            codiAssert(1 == maxOutputArgs);

            // Only called from reuse tape.
            Identifier lhsIdentifier = pointers.lhsIdentifiers[0];
            primalVector[lhsIdentifier] = pointers.oldLhsValues[0];
          }

          /// \copydoc codi::StatementEvaluatorTapeInterface::StatementCallGenerator::evaluate()
          CODI_INLINE static void evaluate(Real* CODI_RESTRICT primalVector, STMT_COMMON_ARGS) {
            Base::internalEvaluate(primalVector, STMT_COMMON_CALL);
          }
      };

      /// IterateInputs implementation
      template<typename Stmt>
      struct StatementCallGenerator<StatementCall::IterateInputs, Stmt>
          : public StatementCallGeneratorBase<Stmt, StatementCallGenerator<StatementCall::IterateInputs, Stmt>> {
        public:
          /// Base class abbreviation.
          using Base = StatementCallGeneratorBase<Stmt, StatementCallGenerator<StatementCall::IterateInputs, Stmt>>;

          /// Callback for identifier iteration.
          using IterCallback = typename LowLevelFunctionEntry<TapeImpl, Real, Identifier>::IterCallback;

          /// \copydoc codi::StatementEvaluatorInnerTapeInterface::StatementCallGenerator::evaluateInner()
          CODI_INLINE static void evaluateInner(size_t const& CODI_RESTRICT maxOutputArgs,
                                                size_t const& CODI_RESTRICT maxActiveArgs,
                                                size_t const& CODI_RESTRICT maxConstantArgs,
                                                size_t& CODI_RESTRICT linearAdjointPos, IterCallback func,
                                                void* userData, STMT_COMMON_ARGS) {
            CODI_UNUSED(linearAdjointPos);

            StatementDataPointers pointers = {};
            pointers.populate(maxOutputArgs, maxActiveArgs, numberOfPassiveArguments, maxConstantArgs, byteData);

            for (Config::ArgumentSize i = 0; i < maxActiveArgs; i += 1) {
              if (pointers.rhsIdentifiers[i] >= (Identifier)Config::MaxArgumentSize) {
                func(&pointers.rhsIdentifiers[i], userData);
              }
            }
          }

          /// \copydoc codi::StatementEvaluatorInnerTapeInterface::StatementCallGenerator::evaluateFull()
          template<typename Func, typename... Args>
          CODI_INLINE static void evaluateFull(Func const& func, Args&&... args) {
            func(std::forward<Args>(args)...);
          }

          /// \copydoc codi::StatementEvaluatorTapeInterface::StatementCallGenerator::evaluate()
          CODI_INLINE static void evaluate(size_t& CODI_RESTRICT linearAdjointPos, IterCallback func, void* userData,
                                           STMT_COMMON_ARGS) {
            Base::internalEvaluate(linearAdjointPos, func, userData, STMT_COMMON_CALL);
          }
      };

      /*******************************************************************************/
      /// IterateOutputs implementation
      template<typename Stmt>
      struct StatementCallGenerator<StatementCall::IterateOutputs, Stmt>
          : public StatementCallGeneratorBase<Stmt, StatementCallGenerator<StatementCall::IterateOutputs, Stmt>> {
        public:
          /// Base class abbreviation.
          using Base = StatementCallGeneratorBase<Stmt, StatementCallGenerator<StatementCall::IterateOutputs, Stmt>>;
          using StaticRhs = typename Base::StaticRhs;  ///< See StatementCallGeneratorBase.

          /// Callback for identifier iteration.
          using IterCallback = typename LowLevelFunctionEntry<TapeImpl, Real, Identifier>::IterCallback;

          /// \copydoc codi::StatementEvaluatorInnerTapeInterface::StatementCallGenerator::evaluateInner()
          CODI_INLINE static void evaluateInner(size_t const& CODI_RESTRICT maxOutputArgs,
                                                size_t const& CODI_RESTRICT maxActiveArgs,
                                                size_t const& CODI_RESTRICT maxConstantArgs,
                                                size_t& CODI_RESTRICT linearAdjointPos, IterCallback func,
                                                void* userData, STMT_COMMON_ARGS) {
            StatementDataPointers pointers = {};
            pointers.populate(maxOutputArgs, maxActiveArgs, numberOfPassiveArguments, maxConstantArgs, byteData);

            for (size_t iLhs = 0; iLhs < maxOutputArgs; iLhs += 1) {
              if (LinearIndexHandling) {
                Identifier lhsIdentifier = linearAdjointPos + 1 + iLhs;
                func(&lhsIdentifier, userData);

                codiAssert(lhsIdentifier == (Identifier)(linearAdjointPos + 1 + iLhs));
              } else {
                func(&pointers.lhsIdentifiers[iLhs], userData);
              }
            }
          }

          /// \copydoc codi::StatementEvaluatorInnerTapeInterface::StatementCallGenerator::evaluateFull()
          template<typename Func, typename... Args>
          CODI_INLINE static void evaluateFull(Func const& func, Args&&... args) {
            func(std::forward<Args>(args)...);
          }

          /// \copydoc codi::StatementEvaluatorTapeInterface::StatementCallGenerator::evaluate()
          CODI_INLINE static void evaluate(size_t& CODI_RESTRICT linearAdjointPos, IterCallback func, void* userData,
                                           STMT_COMMON_ARGS) {
            Base::internalEvaluate(linearAdjointPos, func, userData, STMT_COMMON_CALL);
          }
      };

      /// @}

    private:
      static void evalJacobianReverse(ADJOINT_VECTOR_TYPE* adjointVector, Gradient lhsAdjoint, Real const* const values,
                                      Identifier const* const rhsIdentifiers) {
#if CODI_VariableAdjointInterfaceInPrimalTapes
        CODI_UNUSED(lhsAdjoint);
        bool const lhsZero = adjointVector->isLhsZero();
#else
        bool const lhsZero = RealTraits::isTotalZero(lhsAdjoint);
#endif

        if (CODI_ENABLE_CHECK(Config::SkipZeroAdjointEvaluation, !lhsZero)) {
          for (size_t pos = 0; pos < size; pos += 1) {
            Real const& jacobian = values[pos];
#if CODI_VariableAdjointInterfaceInPrimalTapes
            adjointVector->updateAdjointWithLhs(rhsIdentifiers[pos], jacobian);
#else
            adjointVector[rhsIdentifiers[pos]] += jacobian * lhsAdjoint;
#endif
          }
        }
      }
  };

#undef STMT_COMMON_ARGS
#undef STMT_COMMON_CALL

#define CREATE_EXPRESSION(size)                                                                      \
  TapeTypes::StatementEvaluator::template createHandle<Impl, JacobianStatementGenerator<Impl, size>, \
                                                       AssignStatement<ActiveType<Impl>, JacobianExpression<size>>>()

  /// Expressions for manual statement pushes.
  template<typename TapeTypes, typename Impl>
  const CODI_DD(typename TapeTypes::EvalHandle,
                CODI_ANY) PrimalValueBaseTape<TapeTypes, Impl>::jacobianExpressions[Config::MaxArgumentSize] = {
      CREATE_EXPRESSION(0),   CREATE_EXPRESSION(1),   CREATE_EXPRESSION(2),   CREATE_EXPRESSION(3),
      CREATE_EXPRESSION(4),   CREATE_EXPRESSION(5),   CREATE_EXPRESSION(6),   CREATE_EXPRESSION(7),
      CREATE_EXPRESSION(8),   CREATE_EXPRESSION(9),   CREATE_EXPRESSION(10),  CREATE_EXPRESSION(11),
      CREATE_EXPRESSION(12),  CREATE_EXPRESSION(13),  CREATE_EXPRESSION(14),  CREATE_EXPRESSION(15),
      CREATE_EXPRESSION(16),  CREATE_EXPRESSION(17),  CREATE_EXPRESSION(18),  CREATE_EXPRESSION(19),
      CREATE_EXPRESSION(20),  CREATE_EXPRESSION(21),  CREATE_EXPRESSION(22),  CREATE_EXPRESSION(23),
      CREATE_EXPRESSION(24),  CREATE_EXPRESSION(25),  CREATE_EXPRESSION(26),  CREATE_EXPRESSION(27),
      CREATE_EXPRESSION(28),  CREATE_EXPRESSION(29),  CREATE_EXPRESSION(30),  CREATE_EXPRESSION(31),
      CREATE_EXPRESSION(32),  CREATE_EXPRESSION(33),  CREATE_EXPRESSION(34),  CREATE_EXPRESSION(35),
      CREATE_EXPRESSION(36),  CREATE_EXPRESSION(37),  CREATE_EXPRESSION(38),  CREATE_EXPRESSION(39),
      CREATE_EXPRESSION(40),  CREATE_EXPRESSION(41),  CREATE_EXPRESSION(42),  CREATE_EXPRESSION(43),
      CREATE_EXPRESSION(44),  CREATE_EXPRESSION(45),  CREATE_EXPRESSION(46),  CREATE_EXPRESSION(47),
      CREATE_EXPRESSION(48),  CREATE_EXPRESSION(49),  CREATE_EXPRESSION(50),  CREATE_EXPRESSION(51),
      CREATE_EXPRESSION(52),  CREATE_EXPRESSION(53),  CREATE_EXPRESSION(54),  CREATE_EXPRESSION(55),
      CREATE_EXPRESSION(56),  CREATE_EXPRESSION(57),  CREATE_EXPRESSION(58),  CREATE_EXPRESSION(59),
      CREATE_EXPRESSION(60),  CREATE_EXPRESSION(61),  CREATE_EXPRESSION(62),  CREATE_EXPRESSION(63),
      CREATE_EXPRESSION(64),  CREATE_EXPRESSION(65),  CREATE_EXPRESSION(66),  CREATE_EXPRESSION(67),
      CREATE_EXPRESSION(68),  CREATE_EXPRESSION(69),  CREATE_EXPRESSION(70),  CREATE_EXPRESSION(71),
      CREATE_EXPRESSION(72),  CREATE_EXPRESSION(73),  CREATE_EXPRESSION(74),  CREATE_EXPRESSION(75),
      CREATE_EXPRESSION(76),  CREATE_EXPRESSION(77),  CREATE_EXPRESSION(78),  CREATE_EXPRESSION(79),
      CREATE_EXPRESSION(80),  CREATE_EXPRESSION(81),  CREATE_EXPRESSION(82),  CREATE_EXPRESSION(83),
      CREATE_EXPRESSION(84),  CREATE_EXPRESSION(85),  CREATE_EXPRESSION(86),  CREATE_EXPRESSION(87),
      CREATE_EXPRESSION(88),  CREATE_EXPRESSION(89),  CREATE_EXPRESSION(90),  CREATE_EXPRESSION(91),
      CREATE_EXPRESSION(92),  CREATE_EXPRESSION(93),  CREATE_EXPRESSION(94),  CREATE_EXPRESSION(95),
      CREATE_EXPRESSION(96),  CREATE_EXPRESSION(97),  CREATE_EXPRESSION(98),  CREATE_EXPRESSION(99),
      CREATE_EXPRESSION(100), CREATE_EXPRESSION(101), CREATE_EXPRESSION(102), CREATE_EXPRESSION(103),
      CREATE_EXPRESSION(104), CREATE_EXPRESSION(105), CREATE_EXPRESSION(106), CREATE_EXPRESSION(107),
      CREATE_EXPRESSION(108), CREATE_EXPRESSION(109), CREATE_EXPRESSION(110), CREATE_EXPRESSION(111),
      CREATE_EXPRESSION(112), CREATE_EXPRESSION(113), CREATE_EXPRESSION(114), CREATE_EXPRESSION(115),
      CREATE_EXPRESSION(116), CREATE_EXPRESSION(117), CREATE_EXPRESSION(118), CREATE_EXPRESSION(119),
      CREATE_EXPRESSION(120), CREATE_EXPRESSION(121), CREATE_EXPRESSION(122), CREATE_EXPRESSION(123),
      CREATE_EXPRESSION(124), CREATE_EXPRESSION(125), CREATE_EXPRESSION(126), CREATE_EXPRESSION(127),
      CREATE_EXPRESSION(128), CREATE_EXPRESSION(129), CREATE_EXPRESSION(130), CREATE_EXPRESSION(131),
      CREATE_EXPRESSION(132), CREATE_EXPRESSION(133), CREATE_EXPRESSION(134), CREATE_EXPRESSION(135),
      CREATE_EXPRESSION(136), CREATE_EXPRESSION(137), CREATE_EXPRESSION(138), CREATE_EXPRESSION(139),
      CREATE_EXPRESSION(140), CREATE_EXPRESSION(141), CREATE_EXPRESSION(142), CREATE_EXPRESSION(143),
      CREATE_EXPRESSION(144), CREATE_EXPRESSION(145), CREATE_EXPRESSION(146), CREATE_EXPRESSION(147),
      CREATE_EXPRESSION(148), CREATE_EXPRESSION(149), CREATE_EXPRESSION(150), CREATE_EXPRESSION(151),
      CREATE_EXPRESSION(152), CREATE_EXPRESSION(153), CREATE_EXPRESSION(154), CREATE_EXPRESSION(155),
      CREATE_EXPRESSION(156), CREATE_EXPRESSION(157), CREATE_EXPRESSION(158), CREATE_EXPRESSION(159),
      CREATE_EXPRESSION(160), CREATE_EXPRESSION(161), CREATE_EXPRESSION(162), CREATE_EXPRESSION(163),
      CREATE_EXPRESSION(164), CREATE_EXPRESSION(165), CREATE_EXPRESSION(166), CREATE_EXPRESSION(167),
      CREATE_EXPRESSION(168), CREATE_EXPRESSION(169), CREATE_EXPRESSION(170), CREATE_EXPRESSION(171),
      CREATE_EXPRESSION(172), CREATE_EXPRESSION(173), CREATE_EXPRESSION(174), CREATE_EXPRESSION(175),
      CREATE_EXPRESSION(176), CREATE_EXPRESSION(177), CREATE_EXPRESSION(178), CREATE_EXPRESSION(179),
      CREATE_EXPRESSION(180), CREATE_EXPRESSION(181), CREATE_EXPRESSION(182), CREATE_EXPRESSION(183),
      CREATE_EXPRESSION(184), CREATE_EXPRESSION(185), CREATE_EXPRESSION(186), CREATE_EXPRESSION(187),
      CREATE_EXPRESSION(188), CREATE_EXPRESSION(189), CREATE_EXPRESSION(190), CREATE_EXPRESSION(191),
      CREATE_EXPRESSION(192), CREATE_EXPRESSION(193), CREATE_EXPRESSION(194), CREATE_EXPRESSION(195),
      CREATE_EXPRESSION(196), CREATE_EXPRESSION(197), CREATE_EXPRESSION(198), CREATE_EXPRESSION(199),
      CREATE_EXPRESSION(200), CREATE_EXPRESSION(201), CREATE_EXPRESSION(202), CREATE_EXPRESSION(203),
      CREATE_EXPRESSION(204), CREATE_EXPRESSION(205), CREATE_EXPRESSION(206), CREATE_EXPRESSION(207),
      CREATE_EXPRESSION(208), CREATE_EXPRESSION(209), CREATE_EXPRESSION(210), CREATE_EXPRESSION(211),
      CREATE_EXPRESSION(212), CREATE_EXPRESSION(213), CREATE_EXPRESSION(214), CREATE_EXPRESSION(215),
      CREATE_EXPRESSION(216), CREATE_EXPRESSION(217), CREATE_EXPRESSION(218), CREATE_EXPRESSION(219),
      CREATE_EXPRESSION(220), CREATE_EXPRESSION(221), CREATE_EXPRESSION(222), CREATE_EXPRESSION(223),
      CREATE_EXPRESSION(224), CREATE_EXPRESSION(225), CREATE_EXPRESSION(226), CREATE_EXPRESSION(227),
      CREATE_EXPRESSION(228), CREATE_EXPRESSION(229), CREATE_EXPRESSION(230), CREATE_EXPRESSION(231),
      CREATE_EXPRESSION(232), CREATE_EXPRESSION(233), CREATE_EXPRESSION(234), CREATE_EXPRESSION(235),
      CREATE_EXPRESSION(236), CREATE_EXPRESSION(237), CREATE_EXPRESSION(238), CREATE_EXPRESSION(239),
      CREATE_EXPRESSION(240), CREATE_EXPRESSION(241), CREATE_EXPRESSION(242), CREATE_EXPRESSION(243),
      CREATE_EXPRESSION(244), CREATE_EXPRESSION(245), CREATE_EXPRESSION(246), CREATE_EXPRESSION(247),
      CREATE_EXPRESSION(248), CREATE_EXPRESSION(249), CREATE_EXPRESSION(250), CREATE_EXPRESSION(251),
      CREATE_EXPRESSION(252)};

#undef CREATE_EXPRESSION
}
