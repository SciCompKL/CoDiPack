#pragma once

#include <algorithm>
#include <cmath>
#include <functional>
#include <type_traits>
#include <utility>

#include "../aux/macros.hpp"
#include "../aux/memberStore.hpp"
#include "../config.h"
#include "../expressions/lhsExpressionInterface.hpp"
#include "../expressions/logic/compileTimeTraversalLogic.hpp"
#include "../expressions/logic/constructStaticContext.hpp"
#include "../expressions/logic/helpers/forEachLeafLogic.hpp"
#include "../expressions/logic/helpers/jacobianComputationLogic.hpp"
#include "../expressions/logic/traversalLogic.hpp"
#include "../traits/expressionTraits.hpp"
#include "aux/primalAdjointVectorAccess.hpp"
#include "commonTapeImplementation.hpp"
#include "data/chunk.hpp"
#include "data/chunkedData.hpp"
#include "indices/indexManagerInterface.hpp"
#include "statementEvaluators/statementEvaluatorInterface.hpp"
#include "statementEvaluators/statementEvaluatorTapeInterface.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Type definitions for the primal value tapes.
   *
   * @tparam _Real                See TapeTypesInterface.
   * @tparam _Gradient            See TapeTypesInterface.
   * @tparam _IndexManager        Index manager for the tape. Needs to implement IndexManagerInterface.
   * @tparam _StatementEvaluator  Statement handle generator. Needs to implement StatementEvaluatorInterface and
   *                              StatementEvaluatorInnerTapeInterface.
   * @tparam _Data                See TapeTypesInterface.
   */
  template<typename _Real, typename _Gradient, typename _IndexManager, template<typename> class _StatementEvaluator,
           template<typename, typename> class _Data>
  struct PrimalValueTapeTypes : public TapeTypesInterface {
    public:

      using Real = CODI_DD(_Real, double);                                              ///< See PrimalValueTapeTypes.
      using Gradient = CODI_DD(_Gradient, double);                                      ///< See PrimalValueTapeTypes.
      using IndexManager = CODI_DD(_IndexManager, CODI_T(IndexManagerInterface<int>));  ///< See PrimalValueTapeTypes.
      using StatementEvaluator = CODI_DD(CODI_T(_StatementEvaluator<Real>),
                                         CODI_T(StatementEvaluatorInterface<double>));  ///< See PrimalValueTapeTypes.
      template<typename Chunk, typename Nested>
      using Data = CODI_DD(CODI_T(_Data<Chunk, Nested>), CODI_T(DataInterface<Nested>));  ///< See PrimalValueTapeTypes.

      using Identifier = typename IndexManager::Index;    ///< See IndexManagerInterface.
      using PassiveReal = RealTraits::PassiveReal<Real>;  ///< Basic computation type.

      constexpr static bool IsLinearIndexHandler = IndexManager::IsLinear;  ///< True if the index manager is linear.
      constexpr static bool IsStaticIndexHandler =
          !IsLinearIndexHandler;  ///< For reuse index management, a static index manager is used.

      using EvalHandle = typename StatementEvaluator::Handle;  ///< Handle type returned by the statement generator.

      /// Statement chunk is either \<argument size, eval handle\> (linear management) or \<lhs identifier,
      /// argument size, overwritten value, eval handle\> (reuse management)
      using StatementChunk =
          typename std::conditional<IsLinearIndexHandler, Chunk2<Config::ArgumentSize, EvalHandle>,
                                    Chunk4<Identifier, Config::ArgumentSize, Real, EvalHandle>>::type;
      using StatementData = Data<StatementChunk, IndexManager>;  ///< Statement data vector.

      using IdentifierChunk = Chunk1<Identifier>;                      ///< Identifiers of statement arguments.
      using RhsIdentifierData = Data<IdentifierChunk, StatementData>;  ///< Rhs identifiers data vector.

      using PassiveValueChunk = Chunk1<Real>;                               ///< Passive values of statement arguments.
      using PassiveValueData = Data<PassiveValueChunk, RhsIdentifierData>;  ///< Passive values data vector.

      using ConstantValueChunk = Chunk1<PassiveReal>;  ///< Constant values of in statement expressions.
      using ConstantValueData = Data<ConstantValueChunk, PassiveValueData>;  ///< Constant values data vector.

      using NestedData = ConstantValueData;  ///< See TapeTypesInterface.
  };

  /**
   * @brief Base class for all standard Primal value tape implementations.
   *
   * This class provides nearly a full implementation of the FullTapeInterface. There are just a few internal methods
   * left which need to be implemented by the final classes. These methods depend significantly on the index management
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
   * @tparam _TapeTypes needs to implement PrimalValueTapeTypes.
   * @tparam _Impl Type of the final implementations
   */
  template<typename _TapeTypes, typename _Impl>
  struct PrimalValueBaseTape : public CommonTapeImplementation<_TapeTypes, _Impl>,
                               public StatementEvaluatorTapeInterface<typename _TapeTypes::Real>,
                               public StatementEvaluatorInnerTapeInterface<typename _TapeTypes::Real> {
    public:

      /// See PrimalValueBaseTape.
      using TapeTypes = CODI_DD(_TapeTypes,
                                CODI_T(PrimalValueTapeTypes<double, double, IndexManagerInterface<int>,
                                                            StatementEvaluatorInterface, DefaultChunkedData>));
      /// See PrimalValueBaseTape.
      using Impl = CODI_DD(_Impl, CODI_T(FullTapeInterface<double, double, int, EmptyPosition>));

      using Base = CommonTapeImplementation<TapeTypes, Impl>;  ///< Base class abbreviation.
      friend Base;  ///< Allow the base class to call protected and private methods.

      using Real = typename TapeTypes::Real;                              ///< See TapeTypesInterface.
      using Gradient = typename TapeTypes::Gradient;                      ///< See TapeTypesInterface.
      using IndexManager = typename TapeTypes::IndexManager;              ///< See TapeTypesInterface.
      using StatementEvaluator = typename TapeTypes::StatementEvaluator;  ///< See PrimalValueTapeTypes.
      using Identifier = typename TapeTypes::Identifier;                  ///< See PrimalValueTapeTypes.

      using EvalHandle = typename TapeTypes::EvalHandle;  ///< See PrimalValueTapeTypes.

      using StatementData = typename TapeTypes::StatementData;          ///< See PrimalValueTapeTypes.
      using RhsIdentifierData = typename TapeTypes::RhsIdentifierData;  ///< See PrimalValueTapeTypes.
      using PassiveValueData = typename TapeTypes::PassiveValueData;    ///< See PrimalValueTapeTypes.
      using ConstantValueData = typename TapeTypes::ConstantValueData;  ///< See PrimalValueTapeTypes.

      using PassiveReal = RealTraits::PassiveReal<Real>;  ///< Basic computation type.

      using NestedPosition = typename ConstantValueData::Position;  ///< See PrimalValueTapeTypes.
      using Position = typename Base::Position;                     ///< See TapeTypesInterface.

      template<typename Adjoint>
      using VectorAccess = PrimalAdjointVectorAccess<Real, Identifier, Adjoint>;  ///< Vector access type generated by this tape.

      static bool constexpr AllowJacobianOptimization = false;  ///< See InternalStatementRecordingTapeInterface.
      static bool constexpr HasPrimalValues = true;             ///< See PrimalEvaluationTapeInterface.
      static bool constexpr LinearIndexHandling =
          TapeTypes::IsLinearIndexHandler;  ///< See IdentifierInformationTapeInterface.
      static bool constexpr RequiresPrimalRestore =
          !TapeTypes::IsLinearIndexHandler;  ///< See PrimalEvaluationTapeInterface.

    protected:

      static EvalHandle const jacobianExpressions[Config::MaxArgumentSize];

      MemberStore<IndexManager, Impl, TapeTypes::IsStaticIndexHandler> indexManager;  ///< Index manager.
      StatementData statementData;          ///< Data stream for statement specific data.
      RhsIdentifierData rhsIdentiferData;   ///< Data stream for argument identifier data.
      PassiveValueData passiveValueData;    ///< Data stream for passive argument value data.
      ConstantValueData constantValueData;  ///< Data stream for constant argument data.

      std::vector<Gradient> adjoints;  ///< Evaluation vector for AD.
      std::vector<Real> primals;       ///< Current state of primal values in the program.
      std::vector<Real> primalsCopy;   ///< Copy of primal values for AD evaluations.

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

      /// Perform a primal evaluation of the tape. Arguments are from the recursive eval methods of the DataInterface.
      template<typename... Args>
      static void internalEvaluatePrimal_Step3_EvalStatements(Args&&... args);

      /// Perform a reverse evaluation of the tape. Arguments are from the recursive eval methods of the DataInterface.
      template<typename... Args>
      static void internalEvaluateReverse_Step3_EvalStatements(Args&&... args);

      /// Reset the primal values to the given position.
      void internalResetPrimalValues(Position const& pos);

      /// Add statement specific data to the data streams.
      void pushStmtData(Identifier const& index, Config::ArgumentSize const& numberOfPassiveArguments,
                        Real const& oldPrimalValue, EvalHandle evalHandle);

      /// @}

    public:

      /// Constructor
      PrimalValueBaseTape()
          : Base(),
            indexManager(Config::MaxArgumentSize),  // reserve first items for passive values
            statementData(Config::ChunkSize),
            rhsIdentiferData(Config::ChunkSize),
            passiveValueData(Config::ChunkSize),
            constantValueData(Config::ChunkSize),
            adjoints(1),  // ensure that adjoint[0] exists, see its use in gradient() const
            primals(0),
            primalsCopy(0) {
        checkPrimalSize(true);

        statementData.setNested(&indexManager.get());
        rhsIdentiferData.setNested(&statementData);
        passiveValueData.setNested(&rhsIdentiferData);
        constantValueData.setNested(&passiveValueData);

        Base::init(&constantValueData);

        Base::options.insert(TapeParameters::AdjointSize);
        Base::options.insert(TapeParameters::ConstantValuesSize);
        Base::options.insert(TapeParameters::LargestIdentifier);
        Base::options.insert(TapeParameters::PassiveValuesSize);
        Base::options.insert(TapeParameters::RhsIdentifiersSize);
        Base::options.insert(TapeParameters::PrimalSize);
        Base::options.insert(TapeParameters::StatementSize);
      }

      /*******************************************************************************/
      /// @name Functions from GradientAccessTapeInterface
      /// @{

      /// \copydoc codi::GradientAccessTapeInterface::gradient(Identifier const&)
      CODI_INLINE Gradient& gradient(Identifier const& identifier) {
        checkAdjointSize(identifier);

        return adjoints[identifier];
      }

      /// \copydoc codi::GradientAccessTapeInterface::gradient(Identifier const&) const
      CODI_INLINE Gradient const& gradient(Identifier const& identifier) const {
        if (identifier > (Identifier)adjoints.size()) {
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

        indexManager.get().freeIndex(identifier);
      }

      /// @}

    protected:

      /// Count all arguments that have non-zero index.
      struct CountActiveArguments : public ForEachLeafLogic<CountActiveArguments> {
        public:

          /// \copydoc codi::ForEachLeafLogic::handleActive
          template<typename Node>
          CODI_INLINE void handleActive(Node const& node, size_t& numberOfActiveArguments) {
            if (CODI_ENABLE_CHECK(Config::CheckZeroIndex, IndexManager::InactiveIndex != node.getIdentifier())) {
              numberOfActiveArguments += 1;
            }
          }
      };

      /// Push all data for each argument.
      struct PushIdentfierPassiveAndConstant : public ForEachLeafLogic<PushIdentfierPassiveAndConstant> {
        public:

          /// \copydoc codi::ForEachLeafLogic::handleActive
          template<typename Node>
          CODI_INLINE void handleActive(Node const& node, RhsIdentifierData& rhsIdentiferData,
                                        PassiveValueData& passiveValueData, ConstantValueData& constantValueData,
                                        size_t& curPassiveArgument) {
            CODI_UNUSED(constantValueData);

            Identifier rhsIndex = node.getIdentifier();
            if (CODI_ENABLE_CHECK(Config::CheckZeroIndex, IndexManager::InactiveIndex == rhsIndex)) {
              rhsIndex = curPassiveArgument;

              curPassiveArgument += 1;
              passiveValueData.pushData(node.getValue());
            }

            rhsIdentiferData.pushData(rhsIndex);
          }

          /// \copydoc codi::ForEachLeafLogic::handleConstant
          template<typename Node>
          CODI_INLINE void handleConstant(Node const& node, RhsIdentifierData& rhsIdentiferData,
                                          PassiveValueData& passiveValueData, ConstantValueData& constantValueData,
                                          size_t& curPassiveArgument) {
            CODI_UNUSED(rhsIdentiferData, passiveValueData, curPassiveArgument);

            constantValueData.pushData(node.getValue());
          }
      };

    public:

      /// @{

      /// \copydoc codi::InternalStatementRecordingTapeInterface::store()
      template<typename Lhs, typename Rhs>
      CODI_INLINE void store(LhsExpressionInterface<Real, Gradient, Impl, Lhs>& lhs,
                             ExpressionInterface<Real, Rhs> const& rhs) {
        if (CODI_ENABLE_CHECK(Config::CheckTapeActivity, cast().isActive())) {
          CountActiveArguments countActiveArguments;
          PushIdentfierPassiveAndConstant pushAll;
          size_t constexpr MaxActiveArgs = ExpressionTraits::NumberOfActiveTypeArguments<Rhs>::value;
          size_t constexpr MaxConstantArgs = ExpressionTraits::NumberOfConstantTypeArguments<Rhs>::value;

          codiAssert(MaxActiveArgs < Config::MaxArgumentSize);
          codiAssert(MaxConstantArgs < Config::MaxArgumentSize);

          size_t activeArguments = 0;
          countActiveArguments.eval(rhs.cast(), activeArguments);

          if (0 != activeArguments) {
            statementData.reserveItems(1);
            rhsIdentiferData.reserveItems(MaxActiveArgs);
            passiveValueData.reserveItems(MaxActiveArgs - activeArguments);
            constantValueData.reserveItems(MaxConstantArgs);

            size_t passiveArguments = 0;
            pushAll.eval(rhs.cast(), rhsIdentiferData, passiveValueData, constantValueData, passiveArguments);

            bool generatedNewIndex = indexManager.get().assignIndex(lhs.cast().getIdentifier());
            checkPrimalSize(generatedNewIndex);

            Real& primalEntry = primals[lhs.cast().getIdentifier()];
            cast().pushStmtData(lhs.cast().getIdentifier(), passiveArguments, primalEntry,
                                StatementEvaluator::template createHandle<Impl, Impl, Rhs>());

            primalEntry = rhs.cast().getValue();
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
      CODI_INLINE void store(LhsExpressionInterface<Real, Gradient, Impl, Lhs>& lhs, PassiveReal const& rhs) {
        indexManager.get().freeIndex(lhs.cast().getIdentifier());

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
        bool generatedNewIndex;
        if (unusedIndex) {
          generatedNewIndex = indexManager.get().assignUnusedIndex(value.cast().getIdentifier());
        } else {
          generatedNewIndex = indexManager.get().assignIndex(value.cast().getIdentifier());
        }
        checkPrimalSize(generatedNewIndex);

        Real& primalEntry = primals[value.cast().getIdentifier()];
        if (TapeTypes::IsLinearIndexHandler) {
          statementData.reserveItems(1);
          cast().pushStmtData(value.cast().getIdentifier(), Config::StatementInputTag, primalEntry,
                              StatementEvaluator::template createHandle<Impl, Impl, Lhs>());
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
      }

      /// \copydoc codi::ReverseTapeInterface::clearAdjoints()
      CODI_INLINE void clearAdjoints() {
        for (Gradient& gradient : adjoints) {
          gradient = Gradient();
        }
      }

      /// \copydoc codi::ReverseTapeInterface::reset()
      CODI_INLINE void reset(bool resetAdjoints = true) {
        for (Real& primal : primals) {
          primal = Real();
        }

        Base::reset(resetAdjoints);
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
        values.addDoubleEntry("Memory allocated", memoryAdjoints, true, true);

        values.addSection("Primal vector");
        values.addUnsignedLongEntry("Number of primals", nPrimals);
        values.addDoubleEntry("Memory allocated", memoryPrimals, true, true);

        values.addSection("Index manager");
        indexManager.get().addToTapeValues(values);

        values.addSection("Statement entries");
        statementData.addToTapeValues(values);
        values.addSection("Rhs identifiers entries");
        rhsIdentiferData.addToTapeValues(values);
        values.addSection("Passive value entries");
        passiveValueData.addToTapeValues(values);
        values.addSection("Constant value entries");
        constantValueData.addToTapeValues(values);

        return values;
      }

      /******************************************************************************
       * Protected helper function for CustomAdjointVectorEvaluationTapeInterface
       */

      /// Select the configured adjoint vector, see codi::Config::VariableAdjointInterfaceInPrimalTapes.
      template<typename Adjoint>
      ADJOINT_VECTOR_TYPE* selectAdjointVector(VectorAccess<Adjoint>* vectorAccess, Adjoint* data) {
        CODI_UNUSED(vectorAccess, data);

#if CODI_VariableAdjointInterfaceInPrimalTapes
        return vectorAccess;
#else
        static_assert(std::is_same<Adjoint, Gradient>::value,
                      "Please enable 'CODI_VariableAdjointInterfacePrimalInPrimalTapes' in order"
                      " to use custom adjoint vectors in the primal value tapes.");

        return data;
#endif
      }

      /// Perform the adjoint update based on the configuration in codi::Config::VariableAdjointInterfaceInPrimalTapes.
      struct IncrementForwardLogic : public JacobianComputationLogic<Real, IncrementForwardLogic> {
        public:

          /// \copydoc JacobianComputationLogic::handleJacobianOnActive()
          template<typename Node>
          CODI_INLINE void handleJacobianOnActive(Node const& node, Real jacobian, Gradient& lhsTangent,
                                                  ADJOINT_VECTOR_TYPE* adjointVector) {
            CODI_UNUSED(lhsTangent);

            if (CODI_ENABLE_CHECK(Config::IgnoreInvalidJacobies, RealTraits::isTotalFinite(jacobian))) {
#if CODI_VariableAdjointInterfaceInPrimalTapes
              adjointVector->updateTangentWithLhs(node.getIdentifier(), jacobian);
#else
              lhsTangent += jacobian * adjointVector[node.getIdentifier()];
#endif
            }
          }
      };

      /// Additional wrapper that triggers compiler optimizations.
      CODI_WRAP_FUNCTION(Wrap_internalEvaluateForward_Step3_EvalStatements,
                         Impl::internalEvaluateForward_Step3_EvalStatements);

      /// Forward evaluation of an inner tape part between two external functions.
      CODI_INLINE static void internalEvaluateForward_Step2_DataExtraction(NestedPosition const& start,
                                                                           NestedPosition const& end, Real* primalData,
                                                                           ADJOINT_VECTOR_TYPE* data,
                                                                           ConstantValueData& constantValueData) {
        Wrap_internalEvaluateForward_Step3_EvalStatements evalFunc{};
        constantValueData.evaluateForward(start, end, evalFunc, primalData, data);
      }

      /// Internal method for the forward evaluation of the whole tape.
      template<bool copyPrimal, typename Adjoint>
      CODI_NO_INLINE void internalEvaluateForward(Position const& start, Position const& end, Adjoint* data) {
        std::vector<Real> primalsCopy(0);
        Real* primalData = primals.data();

        if (copyPrimal) {
          primalsCopy = primals;
          primalData = primalsCopy.data();
        }

        VectorAccess<Adjoint> vectorAccess(data, primalData);

        ADJOINT_VECTOR_TYPE* dataVector = selectAdjointVector(&vectorAccess, data);

        Base::internalEvaluateForward_Step1_ExtFunc(start, end, internalEvaluateForward_Step2_DataExtraction,
                                                    &vectorAccess, primalData, dataVector, constantValueData);
      }

      /// Perform the adjoint update based on the configuration in codi::Config::VariableAdjointInterfaceInPrimalTapes.
      struct IncrementReversalLogic : public JacobianComputationLogic<Real, IncrementReversalLogic> {
        public:

          /// See IncrementReversalLogic.
          template<typename Node>
          CODI_INLINE void handleJacobianOnActive(Node const& node, Real jacobian, Gradient const& lhsAdjoint,
                                                  ADJOINT_VECTOR_TYPE* adjointVector) {
            CODI_UNUSED(lhsAdjoint);

            if (CODI_ENABLE_CHECK(Config::IgnoreInvalidJacobies, RealTraits::isTotalFinite(jacobian))) {
#if CODI_VariableAdjointInterfaceInPrimalTapes
              adjointVector->updateAdjointWithLhs(node.getIdentifier(), jacobian);
#else
              adjointVector[node.getIdentifier()] += jacobian * lhsAdjoint;
#endif
            }
          }
      };

      /// Additional wrapper that triggers compiler optimizations.
      CODI_WRAP_FUNCTION(Wrap_internalEvaluateReverse_Step3_EvalStatements,
                         Impl::internalEvaluateReverse_Step3_EvalStatements);

      /// Reverse evaluation of an inner tape part between two external functions.
      CODI_INLINE static void internalEvaluateReverse_Step2_DataExtraction(NestedPosition const& start,
                                                                           NestedPosition const& end, Real* primalData,
                                                                           ADJOINT_VECTOR_TYPE* data,
                                                                           ConstantValueData& constantValueData) {
        Wrap_internalEvaluateReverse_Step3_EvalStatements evalFunc;
        constantValueData.evaluateReverse(start, end, evalFunc, primalData, data);
      }

      /// Internal method for the reverse evaluation of the whole tape.
      template<bool copyPrimal, typename Adjoint>
      CODI_INLINE void internalEvaluateReverse(Position const& start, Position const& end, Adjoint* data) {
        Real* primalData = primals.data();

        if (copyPrimal) {
          primalsCopy = primals;
          primalData = primalsCopy.data();
        }

        VectorAccess<Adjoint> vectorAccess(data, primalData);

        ADJOINT_VECTOR_TYPE* dataVector = selectAdjointVector(&vectorAccess, data);

        Base::internalEvaluateReverse_Step1_ExtFunc(start, end, internalEvaluateReverse_Step2_DataExtraction,
                                                    &vectorAccess, primalData, dataVector, constantValueData);
      }

    public:

      /// @name Functions from CustomAdjointVectorEvaluationTapeInterface
      /// @{

      using Base::evaluate;

      /// \copydoc codi::CustomAdjointVectorEvaluationTapeInterface::evaluate()
      template<typename Adjoint>
      CODI_INLINE void evaluate(Position const& start, Position const& end, Adjoint* data) {
        internalEvaluateReverse<!TapeTypes::IsLinearIndexHandler>(start, end, data);
      }

      /// \copydoc codi::CustomAdjointVectorEvaluationTapeInterface::evaluateForward()
      template<typename Adjoint>
      CODI_INLINE void evaluateForward(Position const& start, Position const& end, Adjoint* data) {
        internalEvaluateForward<!TapeTypes::IsLinearIndexHandler>(start, end, data);
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
      }

      /// \copydoc codi::DataManagementTapeInterface::deleteAdjointVector()
      void deleteAdjointVector() {
        adjoints.resize(1);
      }

      /// \copydoc codi::DataManagementTapeInterface::getParameter()
      size_t getParameter(TapeParameters parameter) const {
        switch (parameter) {
          case TapeParameters::AdjointSize:
            return adjoints.size();
            break;
          case TapeParameters::ConstantValuesSize:
            return constantValueData.getDataSize();
            break;
          case TapeParameters::LargestIdentifier:
            return indexManager.get().getLargestCreatedIndex();
            break;
          case TapeParameters::PassiveValuesSize:
            return passiveValueData.getDataSize();
            break;
          case TapeParameters::RhsIdentifiersSize:
            return rhsIdentiferData.getDataSize();
            break;
          case TapeParameters::PrimalSize:
            return primals.size();
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
          case TapeParameters::ConstantValuesSize:
            constantValueData.resize(value);
            break;
          case TapeParameters::LargestIdentifier:
            CODI_EXCEPTION("Tried to set a get only option.");
            break;
          case TapeParameters::PassiveValuesSize:
            passiveValueData.resize(value);
            break;
          case TapeParameters::RhsIdentifiersSize:
            rhsIdentiferData.resize(value);
            break;
          case TapeParameters::PrimalSize:
            primals.resize(value);
            break;
          case TapeParameters::StatementSize:
            return statementData.resize(value);
          default:
            Base::setParameter(parameter, value);
            break;
        }
      }

      /// \copydoc codi::DataManagementTapeInterface::createVectorAccess()
      VectorAccess<Gradient>* createVectorAccess() {
        return createCustomVectorAccess(adjoints.data());
      }

      /// \copydoc codi::DataManagementTapeInterface::createCustomVectorAccess()
      template<typename Adjoint>
      VectorAccess<Adjoint>* createCustomVectorAccess(Adjoint* data) {
        return new VectorAccess<Adjoint>(data, primals.data());
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
        return internalRegisterInput(value, true);
      }

      /// @}
      /*******************************************************************************/
      /// @name Functions from ForwardEvaluationTapeInterface
      /// @{

      using Base::evaluateForward;

      /// \copydoc codi::ForwardEvaluationTapeInterface::evaluateForward()
      void evaluateForward(Position const& start, Position const& end) {
        checkAdjointSize(indexManager.get().getLargestCreatedIndex());

        cast().evaluateForward(start, end, adjoints.data());
      }

      /// @}
      /******************************************************************************
       * Protected helper function for ManualStatementPushTapeInterface
       */

    protected:

      /// Implements StatementEvaluatorTapeInterface and StatementEvaluatorInnerTapeInterface
      /// @tparam _size  Number of arguments.
      template<size_t _size>
      struct JacobianStatementGenerator {
        public:

          static size_t constexpr size = _size;  ///< See JacobianStatementGenerator

          /*******************************************************************************/
          /// @name Implementation of StatementEvaluatorTapeInterface
          /// @{

          /// Throws exception.
          template<typename Expr, typename... Args>
          static Real statementEvaluateForward(Args&&... args) {
            CODI_UNUSED(args...);
            CODI_EXCEPTION("Forward evaluation of jacobian statement not possible.");
          }

          /// Throws exception.
          template<typename Expr, typename... Args>
          static Real statementEvaluatePrimal(Args&&... args) {
            CODI_UNUSED(args...);
            CODI_EXCEPTION("Primal evaluation of jacobian statement not possible.");
          }

          /// \copydoc codi::StatementEvaluatorTapeInterface::statementEvaluateReverse
          template<typename Expr>
          static void statementEvaluateReverse(Real* primalVector, ADJOINT_VECTOR_TYPE* adjointVector,
                                               Gradient lhsAdjoint, Config::ArgumentSize numberOfPassiveArguments,
                                               size_t& curConstantPos, PassiveReal const* const constantValues,
                                               size_t& curPassivePos, Real const* const passiveValues,
                                               size_t& curRhsIdentifiersPos, Identifier const* const rhsIdentifiers) {
            CODI_UNUSED(primalVector, curConstantPos, constantValues);

            size_t endPos = curRhsIdentifiersPos - numberOfPassiveArguments;

            bool const lhsZero = evalJacobianReverse(adjointVector, lhsAdjoint, curPassivePos, passiveValues,
                                                     curRhsIdentifiersPos, rhsIdentifiers, endPos);

            if (Config::SkipZeroAdjointEvaluation && lhsZero) {
              curPassivePos -= numberOfPassiveArguments;
              curRhsIdentifiersPos -= numberOfPassiveArguments;
            }
          }

          /// @}
          /*******************************************************************************/
          /// @name Implementation of StatementEvaluatorInnerTapeInterface
          /// @{

          /// Throws exception.
          template<typename Expr, typename... Args>
          static Real statementEvaluateForwardInner(Args&&... args) {
            CODI_UNUSED(args...);
            CODI_EXCEPTION("Forward evaluation of jacobian statement not possible.");

            return Real();
          }

          /// Throws exception.
          template<typename Expr, typename... Args>
          static Real statementEvaluatePrimalInner(Args&&... args) {
            CODI_UNUSED(args...);
            CODI_EXCEPTION("Primal evaluation of jacobian statement not possible.");

            return Real();
          }

          /// \copydoc codi::StatementEvaluatorInnerTapeInterface::statementEvaluateReverseInner
          template<typename Expr>
          static void statementEvaluateReverseInner(Real* primalVector, ADJOINT_VECTOR_TYPE* adjointVector,
                                                    Gradient lhsAdjoint, size_t& curConstantPos,
                                                    PassiveReal const* const constantValues,
                                                    size_t& curRhsIdentifiersPos,
                                                    Identifier const* const rhsIdentifiers) {
            CODI_UNUSED(primalVector, curConstantPos, constantValues);

            size_t passivePos = size;
            size_t rhsPos = curRhsIdentifiersPos + size;
            size_t endPos = curRhsIdentifiersPos;

            evalJacobianReverse(adjointVector, lhsAdjoint, passivePos, primalVector, rhsPos, rhsIdentifiers, endPos);
          }

          /// @}

        private:

          static bool evalJacobianReverse(ADJOINT_VECTOR_TYPE* adjointVector, Gradient lhsAdjoint,
                                          size_t& curPassivePos, Real const* const passiveValues,
                                          size_t& curRhsIdentifiersPos, Identifier const* const rhsIdentifiers,
                                          size_t endRhsIdentifiersPos) {
#if CODI_VariableAdjointInterfaceInPrimalTapes
            bool const lhsZero = adjointVector->isLhsZero();
#else
            bool const lhsZero = RealTraits::isTotalZero(lhsAdjoint);
#endif

            if (CODI_ENABLE_CHECK(Config::SkipZeroAdjointEvaluation, !lhsZero)) {
              while (curRhsIdentifiersPos > endRhsIdentifiersPos) {
                curPassivePos -= 1;
                curRhsIdentifiersPos -= 1;

                Real const& jacobian = passiveValues[curPassivePos];
#if CODI_VariableAdjointInterfaceInPrimalTapes
                adjointVector->updateAdjointWithLhs(rhsIdentifiers[curRhsIdentifiersPos], jacobian);
#else
                adjointVector[rhsIdentifiers[curRhsIdentifiersPos]] += jacobian * lhsAdjoint;
#endif
              }
            }

            return lhsZero;
          }
      };

    public:

      /*******************************************************************************/
      /// @name Functions from ManualStatementPushTapeInterface
      /// @{

      /// \copydoc codi::ManualStatementPushTapeInterface::pushJacobiManual()
      void pushJacobiManual(Real const& jacobi, Real const& value, Identifier const& index) {
        CODI_UNUSED(value);

        passiveValueData.pushData(jacobi);
        rhsIdentiferData.pushData(index);
      }

      /// \copydoc codi::ManualStatementPushTapeInterface::storeManual()
      void storeManual(Real const& lhsValue, Identifier& lhsIndex, Config::ArgumentSize const& size) {
        CODI_UNUSED(lhsValue);

        codiAssert(size < Config::MaxArgumentSize);

        statementData.reserveItems(1);
        rhsIdentiferData.reserveItems(size);
        passiveValueData.reserveItems(size);

        indexManager.get().assignIndex(lhsIndex);
        Real& primalEntry = primals[lhsIndex];
        cast().pushStmtData(lhsIndex, size, primalEntry, PrimalValueBaseTape::jacobianExpressions[size]);

        primalEntry = lhsValue;
      }

      /// @}
      /*******************************************************************************/
      /// @name Functions from PositionalEvaluationTapeInterface
      /// @{

      /// \copydoc codi::PositionalEvaluationTapeInterface::evaluate()
      CODI_INLINE void evaluate(Position const& start, Position const& end) {
        checkAdjointSize(indexManager.get().getLargestCreatedIndex());

        evaluate(start, end, adjoints.data());
      }

      /// \copydoc codi::PositionalEvaluationTapeInterface::resetTo()
      CODI_INLINE void resetTo(Position const& pos) {
        cast().internalResetPrimalValues(pos);

        Base::resetTo(pos);
      }

      /// @}
      /*******************************************************************************/
      /// @name Functions from PreaccumulationEvaluationTapeInterface
      /// @{

      /// \copydoc codi::PreaccumulationEvaluationTapeInterface::evaluateKeepState()
      void evaluateKeepState(Position const& start, Position const& end) {
        checkAdjointSize(indexManager.get().getLargestCreatedIndex());

        internalEvaluateReverse<false>(start, end, adjoints.data());

        if (!TapeTypes::IsLinearIndexHandler) {
          evaluatePrimal(end, start);
        }
      }

      /// \copydoc codi::PreaccumulationEvaluationTapeInterface::evaluateForwardKeepState()
      void evaluateForwardKeepState(Position const& start, Position const& end) {
        checkAdjointSize(indexManager.get().getLargestCreatedIndex());

        if (!TapeTypes::IsLinearIndexHandler) {
          cast().internalResetPrimalValues(end);
        }

        internalEvaluateForward<false>(start, end, adjoints.data());
      }

    protected:

      /******************************************************************************
       * Protected helper function for PrimalEvaluationTapeInterface
       */

      /// Wrapper helper for improved compiler optimizations.
      CODI_WRAP_FUNCTION(Wrap_internalEvaluatePrimal_Step3_EvalStatements,
                         Impl::internalEvaluatePrimal_Step3_EvalStatements);

      /// Start for primal evaluation between external function.
      CODI_INLINE static void internalEvaluatePrimal_Step2_DataExtraction(NestedPosition const& start,
                                                                          NestedPosition const& end, Real* primalData,
                                                                          ConstantValueData& constantValueData) {
        Wrap_internalEvaluatePrimal_Step3_EvalStatements evalFunc{};
        constantValueData.evaluateForward(start, end, evalFunc, primalData);
      }

    public:

      /// @}
      /*******************************************************************************/
      /// @name Functions from PrimalEvaluationTapeInterface
      /// @{

      using Base::evaluatePrimal;

      /// \copydoc codi::PrimalEvaluationTapeInterface::evaluatePrimal()
      CODI_NO_INLINE void evaluatePrimal(Position const& start, Position const& end) {
        // TODO: implement primal value only accessor
        PrimalAdjointVectorAccess<Real, Identifier, Gradient> primalAdjointAccess(adjoints.data(), primals.data());

        Base::internalEvaluatePrimal_Step1_ExtFunc(start, end,
                                                   PrimalValueBaseTape::internalEvaluatePrimal_Step2_DataExtraction,
                                                   &primalAdjointAccess, primals.data(), constantValueData);
      }

      /// \copydoc codi::PrimalEvaluationTapeInterface::primal(Identifier const&)
      Real& primal(Identifier const& identifier) {
        return primals[identifier];
      }

      /// \copydoc codi::PrimalEvaluationTapeInterface::primal(Identifier const&) const
      Real const& primal(Identifier const& identifier) const {
        return primals[identifier];
      }

      /// @}
      /*******************************************************************************/
      /// @name Function from StatementEvaluatorInnerTapeInterface
      /// @{

      /// \copydoc codi::StatementEvaluatorInnerTapeInterface::statementEvaluateForwardInner()
      template<typename Rhs>
      static Real statementEvaluateForwardInner(Real* primalVector, ADJOINT_VECTOR_TYPE* adjointVector,
                                                Gradient& lhsTangent, size_t& curConstantPos,
                                                PassiveReal const* const constantValues, size_t& curRhsIdentifiersPos,
                                                Identifier const* const rhsIdentifiers) {
        using Constructor = ConstructStaticContextLogic<Rhs, Impl, 0, 0>;
        using StaticRhs = typename Constructor::ResultType;

        StaticRhs staticsRhs = Constructor::construct(primalVector, &rhsIdentifiers[curRhsIdentifiersPos],
                                                      &constantValues[curConstantPos]);

        IncrementForwardLogic incrementForward;

        incrementForward.eval(staticsRhs, Real(1.0), lhsTangent, adjointVector);
        return staticsRhs.getValue();
      }

      /// \copydoc codi::StatementEvaluatorInnerTapeInterface::statementEvaluateForwardFull()
      template<typename Func>
      static Real statementEvaluateForwardFull(Func const& evalInner, size_t const& maxActiveArgs,
                                               size_t const& maxConstantArgs, Real* primalVector,
                                               ADJOINT_VECTOR_TYPE* adjointVector, Gradient& lhsTangent,
                                               Config::ArgumentSize numberOfPassiveArguments, size_t& curConstantPos,
                                               PassiveReal const* const constantValues, size_t& curPassivePos,
                                               Real const* const passiveValues, size_t& curRhsIdentifiersPos,
                                               Identifier const* const rhsIdentifiers) {
        for (Config::ArgumentSize curPos = 0; curPos < numberOfPassiveArguments; curPos += 1) {
          primalVector[curPos] = passiveValues[curPassivePos + curPos];
        }

        Real ret = evalInner(primalVector, adjointVector, lhsTangent, curConstantPos, constantValues,
                             curRhsIdentifiersPos, rhsIdentifiers);

        // adapt vector positions
        curConstantPos += maxConstantArgs;
        curPassivePos += numberOfPassiveArguments;
        curRhsIdentifiersPos += maxActiveArgs;

        return ret;
      }

      /// \copydoc codi::StatementEvaluatorInnerTapeInterface::statementEvaluatePrimalInner()
      template<typename Rhs>
      static Real statementEvaluatePrimalInner(Real* primalVector, size_t& curConstantPos,
                                               PassiveReal const* const constantValues, size_t& curRhsIdentifiersPos,
                                               Identifier const* const rhsIdentifiers) {
        using Constructor = ConstructStaticContextLogic<Rhs, Impl, 0, 0>;
        using StaticRhs = typename Constructor::ResultType;

        StaticRhs staticsRhs = Constructor::construct(primalVector, &rhsIdentifiers[curRhsIdentifiersPos],
                                                      &constantValues[curConstantPos]);

        return staticsRhs.getValue();
      }

      /// \copydoc codi::StatementEvaluatorInnerTapeInterface::statementEvaluatePrimalFull()
      template<typename Func>
      static Real statementEvaluatePrimalFull(Func const& evalInner, size_t const& maxActiveArgs,
                                              size_t const& maxConstantArgs, Real* primalVector,
                                              Config::ArgumentSize numberOfPassiveArguments, size_t& curConstantPos,
                                              PassiveReal const* const constantValues, size_t& curPassivePos,
                                              Real const* const passiveValues, size_t& curRhsIdentifiersPos,
                                              Identifier const* const rhsIdentifiers) {
        for (Config::ArgumentSize curPos = 0; curPos < numberOfPassiveArguments; curPos += 1) {
          primalVector[curPos] = passiveValues[curPassivePos + curPos];
        }

        Real ret = evalInner(primalVector, curConstantPos, constantValues, curRhsIdentifiersPos, rhsIdentifiers);

        // adapt vector positions
        curConstantPos += maxConstantArgs;
        curPassivePos += numberOfPassiveArguments;
        curRhsIdentifiersPos += maxActiveArgs;

        return ret;
      }

      /// \copydoc codi::StatementEvaluatorInnerTapeInterface::statementEvaluateReverseInner()
      template<typename Rhs>
      CODI_INLINE static void statementEvaluateReverseInner(Real* primalVector, ADJOINT_VECTOR_TYPE* adjointVector,
                                                            Gradient lhsAdjoint, size_t& curConstantPos,
                                                            PassiveReal const* const constantValues,
                                                            size_t& curRhsIdentifiersPos,
                                                            Identifier const* const rhsIdentifiers) {
        using Constructor = ConstructStaticContextLogic<Rhs, Impl, 0, 0>;
        using StaticRhs = typename Constructor::ResultType;

        StaticRhs staticsRhs = Constructor::construct(primalVector, &rhsIdentifiers[curRhsIdentifiersPos],
                                                      &constantValues[curConstantPos]);

        IncrementReversalLogic incrementReverse;

        incrementReverse.eval(staticsRhs, Real(1.0), const_cast<Gradient const&>(lhsAdjoint), adjointVector);
      }

      /// \copydoc codi::StatementEvaluatorInnerTapeInterface::statementEvaluateReverseFull()
      template<typename Func>
      CODI_INLINE static void statementEvaluateReverseFull(
          Func const& evalInner, size_t const& maxActiveArgs, size_t const& maxConstantArgs, Real* primalVector,
          ADJOINT_VECTOR_TYPE* adjointVector, Gradient lhsAdjoint, Config::ArgumentSize numberOfPassiveArguments,
          size_t& curConstantPos, PassiveReal const* const constantValues, size_t& curPassivePos,
          Real const* const passiveValues, size_t& curRhsIdentifiersPos, Identifier const* const rhsIdentifiers) {
        // Adapt vector positions
        curConstantPos -= maxConstantArgs;
        curPassivePos -= numberOfPassiveArguments;
        curRhsIdentifiersPos -= maxActiveArgs;

        if (CODI_ENABLE_CHECK(Config::SkipZeroAdjointEvaluation, !RealTraits::isTotalZero(lhsAdjoint))) {
          for (Config::ArgumentSize curPos = 0; curPos < numberOfPassiveArguments; curPos += 1) {
            primalVector[curPos] = passiveValues[curPassivePos + curPos];
          }

          evalInner(primalVector, adjointVector, lhsAdjoint, curConstantPos, constantValues, curRhsIdentifiersPos,
                    rhsIdentifiers);
        }
      }

      /// @}
      /*******************************************************************************/
      /// @name Function from StatementEvaluatorTapeInterface
      /// @{

      /// \copydoc codi::StatementEvaluatorTapeInterface::statementEvaluateForward()
      template<typename Rhs>
      static Real statementEvaluateForward(Real* primalVector, ADJOINT_VECTOR_TYPE* adjointVector, Gradient& lhsTangent,
                                           Config::ArgumentSize numberOfPassiveArguments, size_t& curConstantPos,
                                           PassiveReal const* const constantValues, size_t& curPassivePos,
                                           Real const* const passiveValues, size_t& curRhsIdentifiersPos,
                                           Identifier const* const rhsIdentifiers) {
        size_t constexpr MaxActiveArgs = ExpressionTraits::NumberOfActiveTypeArguments<Rhs>::value;
        size_t constexpr MaxConstantArgs = ExpressionTraits::NumberOfConstantTypeArguments<Rhs>::value;

        return statementEvaluateForwardFull(statementEvaluateForwardInner<Rhs>, MaxActiveArgs, MaxConstantArgs,
                                            primalVector, adjointVector, lhsTangent, numberOfPassiveArguments,
                                            curConstantPos, constantValues, curPassivePos, passiveValues,
                                            curRhsIdentifiersPos, rhsIdentifiers);
      }

      /// \copydoc codi::StatementEvaluatorTapeInterface::statementEvaluatePrimal()
      template<typename Rhs>
      static Real statementEvaluatePrimal(Real* primalVector, Config::ArgumentSize numberOfPassiveArguments,
                                          size_t& curConstantPos, PassiveReal const* const constantValues,
                                          size_t& curPassivePos, Real const* const passiveValues,
                                          size_t& curRhsIdentifiersPos, Identifier const* const rhsIdentifiers) {
        size_t constexpr MaxActiveArgs = ExpressionTraits::NumberOfActiveTypeArguments<Rhs>::value;
        size_t constexpr MaxConstantArgs = ExpressionTraits::NumberOfConstantTypeArguments<Rhs>::value;

        return statementEvaluatePrimalFull(statementEvaluatePrimalInner<Rhs>, MaxActiveArgs, MaxConstantArgs,
                                           primalVector, numberOfPassiveArguments, curConstantPos, constantValues,
                                           curPassivePos, passiveValues, curRhsIdentifiersPos, rhsIdentifiers);
      }

      /// \copydoc codi::StatementEvaluatorTapeInterface::statementEvaluateReverse()
      template<typename Rhs>
      CODI_INLINE static void statementEvaluateReverse(Real* primalVector, ADJOINT_VECTOR_TYPE* adjointVector,
                                                       Gradient lhsAdjoint,
                                                       Config::ArgumentSize numberOfPassiveArguments,
                                                       size_t& curConstantPos, PassiveReal const* const constantValues,
                                                       size_t& curPassivePos, Real const* const passiveValues,
                                                       size_t& curRhsIdentifiersPos,
                                                       Identifier const* const rhsIdentifiers) {
        size_t constexpr maxActiveArgs = ExpressionTraits::NumberOfActiveTypeArguments<Rhs>::value;
        size_t constexpr maxConstantArgs = ExpressionTraits::NumberOfConstantTypeArguments<Rhs>::value;
        statementEvaluateReverseFull(statementEvaluateReverseInner<Rhs>, maxActiveArgs, maxConstantArgs, primalVector,
                                     adjointVector, lhsAdjoint, numberOfPassiveArguments, curConstantPos,
                                     constantValues, curPassivePos, passiveValues, curRhsIdentifiersPos,
                                     rhsIdentifiers);
      }

    private:

      CODI_INLINE void checkAdjointSize(Identifier const& identifier) {
        if (identifier >= (Identifier)adjoints.size()) {
          resizeAdjointsVector();
        }
      }

      CODI_INLINE void checkPrimalSize(bool generatedNewIndex) {
        if (TapeTypes::IsLinearIndexHandler) {
          if (indexManager.get().getLargestCreatedIndex() >= (Identifier)primals.size()) {
            resizePrimalVector(primals.size() + Config::ChunkSize);
          }
        } else {
          if (generatedNewIndex) {
            resizePrimalVector(indexManager.get().getLargestCreatedIndex() + 1);
          }
        }
      }

      CODI_NO_INLINE void resizeAdjointsVector() {
        adjoints.resize(indexManager.get().getLargestCreatedIndex() + 1);
      }

      CODI_NO_INLINE void resizePrimalVector(size_t newSize) {
        primals.resize(newSize);
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

#define CREATE_EXPRESSION(size)                                                                                        \
  TapeTypes::StatementEvaluator::template createHandle<Impl, typename Impl::template JacobianStatementGenerator<size>, \
                                                       JacobianExpression<size>>()

  /// Expressions for manual statement pushes.
  template<typename TapeTypes, typename Impl>
  const typename TapeTypes::EvalHandle
      PrimalValueBaseTape<TapeTypes, Impl>::jacobianExpressions[Config::MaxArgumentSize] = {
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
          CREATE_EXPRESSION(252), CREATE_EXPRESSION(253)};

#undef CREATE_EXPRESSION
}
