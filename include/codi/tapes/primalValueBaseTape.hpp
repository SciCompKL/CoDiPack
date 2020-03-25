#pragma once

#include <algorithm>
#include <functional>
#include <type_traits>

#include "../aux/macros.h"
#include "../aux/memberStore.hpp"
#include "../config.h"
#include "../expressions/lhsExpressionInterface.hpp"
#include "../expressions/logic/compileTimeTraversalLogic.hpp"
#include "../expressions/logic/traversalLogic.hpp"
#include "../expressions/logic/constructStaticContext.hpp"
#include "../traits/expressionTraits.hpp"
#include "aux/primalAdjointVectorAccess.hpp"
#include "data/chunk.hpp"
#include "data/chunkVector.hpp"
#include "indices/indexManagerInterface.hpp"
#include "commonTapeImplementation.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  template<typename _Real, typename _Gradient, typename _IndexManager>
  struct PrimalValueTapeTypes : public TapeTypesInterface {
    public:

      using Real = DECLARE_DEFAULT(_Real, double);
      using Gradient = DECLARE_DEFAULT(_Gradient, double);
      using IndexManager = DECLARE_DEFAULT(_IndexManager, TEMPLATE(IndexManagerInterface<int>));

      using Identifier = typename IndexManager::Index;
      using PassiveReal = PassiveRealType<Real>;

      constexpr static bool IsLinearIndexHandler = IndexManager::IsLinear;
      constexpr static bool IsStaticIndexHandler = !IsLinearIndexHandler;

      using EvalPointer = void (*)(
          Real* primalVector, ADJOINT_VECTOR_TYPE* adjointVector,
          Gradient lhsAdjoint, Config::ArgumentSize numberOfPassiveArguments,
          size_t& curConstantPos, PassiveReal const* const constantValues,
          size_t& curPassivePos, Real const* const passiveValues,
          size_t& curRhsIdentifiersPos, Identifier const* const rhsIdentifiers);

      using StatementChunk = typename std::conditional<
                                IsLinearIndexHandler,
                                Chunk2<Config::ArgumentSize, EvalPointer>,
                                Chunk4<Identifier, Config::ArgumentSize, Real, EvalPointer>
                              >::type;
      using StatementVector = ChunkVector<StatementChunk, IndexManager>;

      using IdentifierChunk = Chunk1<Identifier>;
      using RhsIdentifierVector = ChunkVector<IdentifierChunk, StatementVector>;

      using PassiveValueChunk = Chunk1<Real>;
      using PassiveValueVector = ChunkVector<PassiveValueChunk, RhsIdentifierVector>;

      using ConstantValueChunk = Chunk1<PassiveReal>;
      using ConstantValueVector = ChunkVector<ConstantValueChunk, PassiveValueVector>;

      using NestedVector = ConstantValueVector;
  };

  template<typename _TapeTypes, typename _Impl>
  struct PrimalValueBaseTape :  public CommonTapeImplementation<_TapeTypes, _Impl> {
    public:

      using TapeTypes = DECLARE_DEFAULT(_TapeTypes, TEMPLATE(PrimalValueTapeTypes<double, double, IndexManagerInterface<int>));
      using Impl = DECLARE_DEFAULT(_Impl, TEMPLATE(ReverseTapeInterface<double, double, int>));

      using Base = CommonTapeImplementation<TapeTypes, Impl>;
      friend Base;

      using Real = typename TapeTypes::Real;
      using Gradient = typename TapeTypes::Gradient;
      using IndexManager = typename TapeTypes::IndexManager;
      using Identifier = typename TapeTypes::Identifier;

      using EvalPointer = typename TapeTypes::EvalPointer;

      using StatementVector = typename TapeTypes::StatementVector;
      using RhsIdentifierVector = typename TapeTypes::RhsIdentifierVector;
      using PassiveValueVector = typename TapeTypes::PassiveValueVector;
      using ConstantValueVector = typename TapeTypes::ConstantValueVector;

      using PassiveReal = PassiveRealType<Real>;

      using NestedPosition = typename ConstantValueVector::Position;
      using Position = typename Base::Position;

      static bool constexpr AllowJacobianOptimization = false;

    protected:

      MemberStore<IndexManager, Impl, TapeTypes::IsStaticIndexHandler> indexManager;
      StatementVector statementVector;
      RhsIdentifierVector rhsIdentiferVector;
      PassiveValueVector passiveValueVector;
      ConstantValueVector constantValueVector;

      std::vector<Gradient> adjoints;
      std::vector<Real> primals;
      std::vector<Real> primalsCopy;

    public:

      CODI_INLINE Impl const& cast() const {
        return static_cast<Impl const&>(*this);
      }

      CODI_INLINE Impl& cast() {
        return static_cast<Impl&>(*this);
      }

      /*******************************************************************************
       * Section: Methods expected in the child class.
       *
       * Description: TODO
       *
       */
    protected:

      template<typename ... Args>
      static void internalEvaluateForwardStack(Args&& ... args);

      template<typename ... Args>
      static void internalEvaluatePrimal(Args&& ... args);

      template<typename ... Args>
      static void internalEvaluateReverseStack(Args&& ... args);

      void internalResetPrimalValues(Position const& pos);

      void pushStmtData(
          Identifier const& index,
          Config::ArgumentSize const& numberOfPassiveArguments,
          Real const& oldPrimalValue,
          EvalPointer evalFunc);

    public:

      PrimalValueBaseTape() :
        Base(),
        indexManager(Config::MaxArgumentSize),
        statementVector(Config::ChunkSize),
        rhsIdentiferVector(Config::ChunkSize),
        passiveValueVector(Config::ChunkSize),
        constantValueVector(Config::ChunkSize),
        adjoints(1),
        primals(0),
        primalsCopy(0)
      {
        checkPrimalSize(true);

        statementVector.setNested(&indexManager.get());
        rhsIdentiferVector.setNested(&statementVector);
        passiveValueVector.setNested(&rhsIdentiferVector);
        constantValueVector.setNested(&passiveValueVector);

        Base::init(&constantValueVector);

        Base::options.insert(ConfigurationOption::AdjointSize);
        Base::options.insert(ConfigurationOption::ConstantValuesSize);
        Base::options.insert(ConfigurationOption::PassiveValuesSize);
        Base::options.insert(ConfigurationOption::RhsIdentifiersSize);
        Base::options.insert(ConfigurationOption::PrimalSize);
        Base::options.insert(ConfigurationOption::StatementSize);
      }

      /*******************************************************************************
       * Section: Functions from GradientAccessInterface
       *
       */

      CODI_INLINE Gradient& gradient(Identifier const& identifier) {
        checkAdjointSize(identifier);

        return adjoints[identifier];
      }

      CODI_INLINE Gradient const& gradient(Identifier const& identifier) const {
        if(identifier > (Identifier)adjoints.size()) {
          return adjoints[0];
        } else {
          return adjoints[identifier];
        }
      }

      /*******************************************************************************
       * Section: Functions from InternalExpressionTapeInterface
       *
       */

      template<typename Real>
      CODI_INLINE void initIdentifier(Real& value, Identifier& identifier) {
        CODI_UNUSED(value);

        identifier = IndexManager::UnusedIndex;
      }

      template<typename Real>
      CODI_INLINE void destroyIdentifier(Real& value, Identifier& identifier) {
        CODI_UNUSED(value);

        indexManager.get().freeIndex(identifier);
      }

    protected:

      struct CountActiveArguments : public TraversalLogic<CountActiveArguments> {
        public:
          template<typename Node>
          CODI_INLINE enableIfLhsExpression<Node> term(Node const& node, size_t& numberOfActiveArguments) {
            ENABLE_CHECK(Config::CheckZeroIndex, 0 != node.getIdentifier()) {

              numberOfActiveArguments += 1;
            }
          }
          using TraversalLogic<CountActiveArguments>::term;
      };

      struct PushIdentfierPassiveAndConstant : public TraversalLogic<PushIdentfierPassiveAndConstant> {
        public:
          template<typename Node>
          CODI_INLINE enableIfLhsExpression<Node> term(
              Node const& node,
              RhsIdentifierVector& rhsIdentiferVector,
              PassiveValueVector& passiveValueVector,
              ConstantValueVector& constantValueVector,
              size_t& curPassiveArgument) {

            CODI_UNUSED(constantValueVector);

            Identifier rhsIndex = node.getIdentifier();
            ENABLE_CHECK(Config::CheckZeroIndex, 0 == rhsIndex) {
              rhsIndex = curPassiveArgument;

              curPassiveArgument += 1;
              passiveValueVector.pushData(node.getValue());
            }

            rhsIdentiferVector.pushData(rhsIndex);
          }

          template<typename Node>
          CODI_INLINE enableIfConstantExpression<Node> term(
              Node const& node,
              RhsIdentifierVector& rhsIdentiferVector,
              PassiveValueVector& passiveValueVector,
              ConstantValueVector& constantValueVector,
              size_t& curPassiveArgument) {

            CODI_UNUSED(rhsIdentiferVector, passiveValueVector, curPassiveArgument);

            constantValueVector.pushData(node.getValue());
          }

          using TraversalLogic<PushIdentfierPassiveAndConstant>::term;
      };

    public:

      template<typename Lhs, typename Rhs>
      CODI_INLINE void store(LhsExpressionInterface<Real, Gradient, Impl, Lhs>& lhs,
                 ExpressionInterface<Real, Rhs> const& rhs) {

        if(cast().isActive()) {
          CountActiveArguments countActiveArguments;
          PushIdentfierPassiveAndConstant pushAll;
          size_t constexpr MaxActiveArgs = MaxNumberOfActiveTypeArguments<Rhs>::value;
          size_t constexpr MaxConstantArgs = MaxNumberOfConstantArguments<Rhs>::value;

          size_t activeArguments = 0;
          countActiveArguments.eval(rhs.cast(), activeArguments);

          if(0 != activeArguments) {

            statementVector.reserveItems(1);
            rhsIdentiferVector.reserveItems(MaxActiveArgs);
            passiveValueVector.reserveItems(MaxActiveArgs - activeArguments);
            constantValueVector.reserveItems(MaxConstantArgs);

            size_t passiveArguments = 0;
            pushAll.eval(rhs.cast(), rhsIdentiferVector, passiveValueVector, constantValueVector, passiveArguments);

            bool generatedNewIndex = indexManager.get().assignIndex(lhs.cast().getIdentifier());
            checkPrimalSize(generatedNewIndex);

            Real& primalEntry = primals[lhs.cast().getIdentifier()];
            cast().pushStmtData(lhs.cast().getIdentifier(), passiveArguments, primalEntry, statementReversal<Rhs>);

            primalEntry = rhs.cast().getValue();
          } else {
            indexManager.get().freeIndex(lhs.cast().getIdentifier());
          }
        } else {
          indexManager.get().freeIndex(lhs.cast().getIdentifier());
        }

        lhs.cast().value() = rhs.cast().getValue();
      }

      template<typename Lhs, typename Rhs>
      CODI_INLINE void store(LhsExpressionInterface<Real, Gradient, Impl, Lhs>& lhs,
                 LhsExpressionInterface<Real, Gradient, Impl, Rhs> const& rhs) {

        if(cast().isActive()) {
          if(IndexManager::AssignNeedsStatement || !Config::AssignOptimization) {
            store<Lhs, Rhs>(lhs, static_cast<ExpressionInterface<Real, Rhs> const&>(rhs));
          } else {
            indexManager.get().copyIndex(lhs.cast().getIdentifier(), rhs.cast().getIdentifier());
          }
        } else {
          indexManager.get().freeIndex(lhs.cast().getIdentifier());
        }

        lhs.cast().value() = rhs.cast().getValue();
      }

      template<typename Lhs>
      CODI_INLINE void store(LhsExpressionInterface<Real, Gradient, Impl, Lhs>& lhs, PassiveReal const& rhs) {
        indexManager.get().freeIndex(lhs.cast().getIdentifier());

        lhs.cast().value() = rhs;
      }

      /*******************************************************************************
       * Section: Functions from ReverseTapeInterface
       *
       */

    protected:

      template<typename Lhs>
      CODI_INLINE Real internalRegisterInput(LhsExpressionInterface<Real, Gradient, Impl, Lhs>& value, bool unusedIndex) {
        bool generatedNewIndex;
        if(unusedIndex) {
          generatedNewIndex = indexManager.get().assignUnusedIndex(value.cast().getIdentifier());
        } else {
          generatedNewIndex = indexManager.get().assignIndex(value.cast().getIdentifier());
        }
        checkPrimalSize(generatedNewIndex);

        Real& primalEntry = primals[value.cast().getIdentifier()];
        if(TapeTypes::IsLinearIndexHandler) {
          statementVector.reserveItems(1);
          cast().pushStmtData(value.cast().getIdentifier(), Config::StatementInputTag, primalEntry, statementReversal<Lhs>);
        }

        Real oldValue = primalEntry;
        primalEntry = value.cast().value();

        return oldValue;
      }

    public:

      template<typename Lhs>
      CODI_INLINE void registerInput(LhsExpressionInterface<Real, Gradient, Impl, Lhs>& value) {
        internalRegisterInput(value, true);
      }

      void clearAdjoints() {
        for(Gradient& gradient : adjoints) {
          gradient = Gradient();
        }
      }

      void reset(bool resetAdjoints = true) {
        for(Real& primal : primals) {
          primal = Real();
        }

        Base::reset(resetAdjoints);
      }

    protected:

      TapeValues internalGetTapeValues() const {
        std::string name;
        if(TapeTypes::IsLinearIndexHandler) {
          name = "CoDi Tape Statistics ( PrimalValueLinearTape )";
        } else {
          name = "CoDi Tape Statistics ( PrimalValueReuseTape )";
        }
        TapeValues values = TapeValues(name);

        size_t nAdjoints      = indexManager.get().getLargestAssignedIndex();
        double memoryAdjoints = static_cast<double>(nAdjoints) * static_cast<double>(sizeof(Gradient)) * TapeValues::BYTE_TO_MB;

        size_t nPrimals = indexManager.get().getLargestAssignedIndex();
        double memoryPrimals = static_cast<double>(nPrimals) * static_cast<double>(sizeof(Real)) * TapeValues::BYTE_TO_MB;

        values.addSection("Adjoint vector");
        values.addUnsignedLongEntry("Number of adjoints", nAdjoints);
        values.addDoubleEntry("Memory allocated", memoryAdjoints, true, true);

        values.addSection("Primal vector");
        values.addUnsignedLongEntry("Number of primals", nPrimals);
        values.addDoubleEntry("Memory allocated", memoryPrimals, true, true);

        values.addSection("Index manager");
        indexManager.get().addToTapeValues(values);

        values.addSection("Statement entries");
        statementVector.addToTapeValues(values);
        values.addSection("Rhs identifiers entries");
        rhsIdentiferVector.addToTapeValues(values);
        values.addSection("Passive value entries");
        passiveValueVector.addToTapeValues(values);
        values.addSection("Constant value entries");
        constantValueVector.addToTapeValues(values);

        return values;
      }

    public:

      using Base::evaluate;

      /*******************************************************************************
       * Section: Function from CustomVectorEvaluationTapeInterface
       *
       */

    protected:

      template<typename Adjoint>
      ADJOINT_VECTOR_TYPE* wrapAdjointVector(VectorAccessInterface<Real, Identifier>* vectorAccess, Adjoint* data) {
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



      struct IncrementReversalLogic : public TraversalLogic<IncrementReversalLogic> {
        public:
          template<typename Node>
          CODI_INLINE enableIfStaticContextActiveType<Node> term(Node const& node, Real jacobian, Gradient const& lhsAdjoint, ADJOINT_VECTOR_TYPE* adjointVector) {
            CODI_UNUSED(lhsAdjoint);

            using std::isfinite;
            ENABLE_CHECK(Config::IgnoreInvalidJacobies, isfinite(jacobian)) {
#if CODI_VariableAdjointInterfaceInPrimalTapes
              adjointVector->updateAdjointWithLhs(node.getIdentifier(), jacobian);
#else
              adjointVector[node.getIdentifier()] += jacobian * lhsAdjoint;
#endif
            }
          }
          using TraversalLogic<IncrementReversalLogic>::term;

          template<size_t LeafNumber, typename Leaf, typename Root>
          CODI_INLINE void link(Leaf const& leaf, Root const& root, Real jacobian, Gradient const& lhsAdjoint, Gradient* adjointVector) {

            Real curJacobian = root.template getJacobian<LeafNumber>() * jacobian;

            this->toNode(leaf, curJacobian, lhsAdjoint, adjointVector);
          }
      };

      template<typename Rhs>
      static void statementReversal(
          Real* primalVector, ADJOINT_VECTOR_TYPE* adjointVector,
          Gradient lhsAdjoint, Config::ArgumentSize numberOfPassiveArguments,
          size_t& curConstantPos, PassiveReal const* const constantValues,
          size_t& curPassivePos, Real const* const passiveValues,
          size_t& curRhsIdentifiersPos, Identifier const* const rhsIdentifiers) {

        using Construtor = ConstructStaticContextlLogic<Rhs, Impl, 0, 0>;
        using StaticRhs = typename Construtor::ResultType;

        size_t constexpr MaxActiveArgs = MaxNumberOfActiveTypeArguments<Rhs>::value;
        size_t constexpr MaxConstantArgs = MaxNumberOfConstantArguments<Rhs>::value;

        // Adapt vector positions
        curConstantPos -= MaxConstantArgs;
        curPassivePos -= numberOfPassiveArguments;
        curRhsIdentifiersPos -= MaxActiveArgs;

        ENABLE_CHECK(Config::SkipZeroAdjointEvaluation, !isTotalZero(lhsAdjoint)) {
          for(Config::ArgumentSize curPos = 0; curPos < numberOfPassiveArguments; curPos += 1) {
            primalVector[curPos] = passiveValues[curPassivePos + curPos];
          }

          StaticRhs staticsRhs = Construtor::construct(
                primalVector,
                &rhsIdentifiers[curRhsIdentifiersPos],
                &constantValues[curConstantPos]);

          IncrementReversalLogic incrementReverse;

          incrementReverse.eval(staticsRhs, 1.0, const_cast<Gradient const&>(lhsAdjoint), adjointVector);
        }
      }

      WRAP_FUNCTION(Wrap_internalEvaluateReverseStack, Impl::internalEvaluateReverseStack);

      CODI_INLINE static void internalEvaluateReverseVector(NestedPosition const& start, NestedPosition const& end,
                                                   Real* primalData,
                                                   ADJOINT_VECTOR_TYPE* data,
                                                   ConstantValueVector& constantValueVector) {
        Wrap_internalEvaluateReverseStack evalFunc;
        constantValueVector.evaluateReverse(start, end, evalFunc, primalData, data);
      }

      WRAP_FUNCTION(Wrap_internalEvaluateReverseVector, internalEvaluateReverseVector);

      template<bool copyPrimal, typename Adjoint>
      CODI_INLINE void internalEvaluateReverse(Position const& start, Position const& end, Adjoint* data) {

        Real* primalData = primals.data();

        if(copyPrimal) {
          primalsCopy = primals;
          primalData = primalsCopy.data();
        }

        AdjointVectorAccess<Real, Identifier, Adjoint> adjointAccess(data);
        PrimalAdjointVectorAccess<Real, Identifier, Adjoint> primalAdjointAccess(data, primalData);

        VectorAccessInterface<Real, Identifier>* vectorAccess;

        if(TapeTypes::IsLinearIndexHandler) {
          vectorAccess = &adjointAccess;
        } else {
          vectorAccess = &primalAdjointAccess;
        }

        ADJOINT_VECTOR_TYPE* dataVector = wrapAdjointVector(vectorAccess, data);

        if(TapeTypes::IsLinearIndexHandler) {
          Wrap_internalEvaluateReverseVector evalFunc{};
          Base::internalEvaluateExtFunc(start, end, evalFunc, vectorAccess,
                                      primalData, dataVector, constantValueVector);
        } else {
          Base::internalEvaluateExtFunc(start, end, internalEvaluateReverseVector, vectorAccess,
                                        primalData, dataVector, constantValueVector);
        }
      }

    public:

      template<typename Adjoint>
      CODI_INLINE void evaluate(Position const& start, Position const& end, Adjoint* data) {
        internalEvaluateReverse<!TapeTypes::IsLinearIndexHandler>(start, end, data);
      }

    protected:

      struct IncrementForwardLogic : public TraversalLogic<IncrementForwardLogic> {
        public:
          template<typename Node>
          CODI_INLINE enableIfStaticContextActiveType<Node> term(Node const& node, Real jacobian, Gradient& lhsTangent, ADJOINT_VECTOR_TYPE* adjointVector) {
            CODI_UNUSED(lhsTangent);

            using std::isfinite;
            ENABLE_CHECK(Config::IgnoreInvalidJacobies, isfinite(jacobian)) {
#if CODI_VariableAdjointInterfaceInPrimalTapes
              adjointVector->updateTangentWithLhs(node.getIdentifier(), jacobian);
#else
              lhsTangent += jacobian * adjointVector[node.getIdentifier()];
#endif
            }
          }
          using TraversalLogic<IncrementForwardLogic>::term;

          template<size_t LeafNumber, typename Leaf, typename Root>
          CODI_INLINE void link(Leaf const& leaf, Root const& root, Real jacobian, Gradient const& lhsTangent, Gradient* adjointVector) {

            Real curJacobian = root.template getJacobian<LeafNumber>() * jacobian;

            this->toNode(leaf, curJacobian, lhsTangent, adjointVector);
          }
      };

      template<typename Rhs>
      static Real statementForward(
          Real* primalVector, ADJOINT_VECTOR_TYPE* adjointVector,
          Gradient& lhsTangent, Config::ArgumentSize numberOfPassiveArguments,
          size_t& curConstantPos, PassiveReal const* const constantValues,
          size_t& curPassivePos, Real const* const passiveValues,
          size_t& curRhsIdentifiersPos, Identifier const* const rhsIdentifiers) {

        using Construtor = ConstructStaticContextlLogic<Rhs, Impl, 0, 0>;
        using StaticRhs = typename Construtor::ResultType;

        size_t constexpr MaxActiveArgs = MaxNumberOfActiveTypeArguments<Rhs>::value;
        size_t constexpr MaxConstantArgs = MaxNumberOfConstantArguments<Rhs>::value;

        for(Config::ArgumentSize curPos = 0; curPos < numberOfPassiveArguments; curPos += 1) {
          primalVector[curPos] = passiveValues[curPassivePos + curPos];

        }

        StaticRhs staticsRhs = Construtor::construct(
              primalVector,
              &rhsIdentifiers[curRhsIdentifiersPos],
              &constantValues[curConstantPos]);

        IncrementForwardLogic incrementForward;

        incrementForward.eval(staticsRhs, 1.0, lhsTangent, adjointVector);

        // Adapt vector positions
        curConstantPos += MaxConstantArgs;
        curPassivePos += numberOfPassiveArguments;
        curRhsIdentifiersPos += MaxActiveArgs;

        return staticsRhs.getValue();
      }

      WRAP_FUNCTION(Wrap_internalEvaluateForwardStack, Impl::internalEvaluateForwardStack);

      CODI_INLINE static void internalEvaluateForwardVector(NestedPosition const& start, NestedPosition const& end,
                                                     Real* primalData,
                                                     ADJOINT_VECTOR_TYPE* data,
                                                     ConstantValueVector& constantValueVector) {


        Wrap_internalEvaluateForwardStack evalFunc{};
        constantValueVector.evaluateForward(start, end, evalFunc, primalData, data);
      }

      WRAP_FUNCTION(Wrap_internalEvaluateForwardVector, internalEvaluateForwardVector);

      template<bool copyPrimal, typename Adjoint>
      CODI_NO_INLINE void internalEvaluateForward(Position const& start, Position const& end, Adjoint* data) {

        std::vector<Real> primalsCopy(0);
        Real* primalData = primals.data();

        if(copyPrimal) {
          primalsCopy = primals;
          primalData = primalsCopy.data();
        }

        AdjointVectorAccess<Real, Identifier, Adjoint> adjointAccess(data);
        PrimalAdjointVectorAccess<Real, Identifier, Adjoint> primalAdjointAccess(data, primalData);

        VectorAccessInterface<Real, Identifier>* vectorAccess;

        if(TapeTypes::IsLinearIndexHandler) {
          vectorAccess = &adjointAccess;
        } else {
          vectorAccess = &primalAdjointAccess;
        }

        ADJOINT_VECTOR_TYPE* dataVector = wrapAdjointVector(vectorAccess, data);

        if(TapeTypes::IsLinearIndexHandler) {
          Wrap_internalEvaluateForwardVector evalFunc{};
          Base::internalEvaluateExtFuncForward(start, end, evalFunc, vectorAccess,
                                               primalData, dataVector, constantValueVector);
        } else {
          Base::internalEvaluateExtFuncForward(start, end, internalEvaluateForwardVector, vectorAccess,
                                               primalData, dataVector, constantValueVector);
        }

      }
    public:

      template<typename Adjoint>
      CODI_INLINE void evaluateForward(Position const& start, Position const& end, Adjoint* data) {
        internalEvaluateForward<!TapeTypes::IsLinearIndexHandler>(start, end, data);
      }

      /*******************************************************************************
       * Section: Function from DataManagementTapeInterface
       *
       */

      CODI_INLINE void swap(Impl& other) {

        // Index manager does not need to be swapped, it is either static or swapped in with the vector data
        // Vectors are swapped recursively in the base class

        std::swap(adjoints, other.adjoints);
        std::swap(primals, other.primals);

        Base::swap(other);
      }

      void deleteAdjointVector() {
        adjoints.resize(1);
      }

      size_t getOption(ConfigurationOption option) const {
        switch (option) {
          case ConfigurationOption::AdjointSize:
            return adjoints.size();
            break;
          case ConfigurationOption::ConstantValuesSize:
            return constantValueVector.getDataSize();
            break;
          case ConfigurationOption::PassiveValuesSize:
            return passiveValueVector.getDataSize();
            break;
          case ConfigurationOption::RhsIdentifiersSize:
            return rhsIdentiferVector.getDataSize();
            break;
        case ConfigurationOption::PrimalSize:
          return primals.size();
          break;
          case ConfigurationOption::StatementSize:
            return statementVector.getDataSize();
          default:
            return Base::getOption(option);
            break;
        }
      }

      void setOption(ConfigurationOption option, size_t value) {
        switch (option) {
          case ConfigurationOption::AdjointSize:
            return adjoints.resize(value);
            break;
          case ConfigurationOption::ConstantValuesSize:
            return constantValueVector.resize(value);
            break;
          case ConfigurationOption::PassiveValuesSize:
            return passiveValueVector.resize(value);
            break;
          case ConfigurationOption::RhsIdentifiersSize:
            return rhsIdentiferVector.resize(value);
            break;
          case ConfigurationOption::PrimalSize:
            return primals.resize(value);
            break;
          case ConfigurationOption::StatementSize:
            return statementVector.resize(value);
          default:
            return Base::setOption(option, value);
            break;
        }
      }

      /*******************************************************************************
       * Section: Function from ExternalFunctionTapeInterface
       *
       */

      template<typename Lhs>
      Real registerExternalFunctionOutput(LhsExpressionInterface<Real, Gradient, Impl, Lhs>& value) {
        return internalRegisterInput(value, true);
      }

      /*******************************************************************************
       * Section: Function from ForwardEvaluationTapeInterface
       *
       */

      using Base::evaluateForward;
      void evaluateForward(Position const& start, Position const& end) {

        checkAdjointSize(indexManager.get().getLargestAssignedIndex());

        cast().evaluateForward(start, end, adjoints.data());
      }

      /*******************************************************************************
       * Section: Function from ManualStatementPushTapeInterface
       *
       */

    protected:

      static void manualStatementReversal(
          Real* primalVector, ADJOINT_VECTOR_TYPE* adjointVector,
          Gradient lhsAdjoint, Config::ArgumentSize numberOfPassiveArguments,
          size_t& curConstantPos, PassiveReal const* const constantValues,
          size_t& curPassivePos, Real const* const passiveValues,
          size_t& curRhsIdentifiersPos, Identifier const* const rhsIdentifiers) {

        CODI_UNUSED(primalVector, curConstantPos, constantValues);

        size_t endPos = curPassivePos - numberOfPassiveArguments;

        #if CODI_VariableAdjointInterfaceInPrimalTapes
          bool const lhsZero = adjointVector->isLhsZero();
        #else
          bool const lhsZero = isTotalZero(lhsAdjoint);
        #endif

        ENABLE_CHECK(Config::SkipZeroAdjointEvaluation, !lhsZero) {
          while(curPassivePos > endPos) {
            curPassivePos -= 1;
            curRhsIdentifiersPos -= 1;

            Real const& jacobian = passiveValues[curPassivePos];
            #if CODI_VariableAdjointInterfaceInPrimalTapes
              adjointVector->updateAdjointWithLhs(rhsIdentifiers[curRhsIdentifiersPos], jacobian);
            #else
              adjointVector[rhsIdentifiers[curRhsIdentifiersPos]] += jacobian * lhsAdjoint;
            #endif

          }
        } else {
          curPassivePos -= numberOfPassiveArguments;
          curRhsIdentifiersPos -= numberOfPassiveArguments;
        }
      }

    public:

      void pushJacobiManual(Real const& jacobi, Real const& value, Gradient const& index) {
        CODI_UNUSED(value);

        passiveValueVector.pushData(jacobi);
        rhsIdentiferVector.pushData(index);
      }

      void storeManual(Real const& lhsValue, Gradient& lhsIndex, Config::ArgumentSize const& size) {
        CODI_UNUSED(lhsValue);

        statementVector.reserveItems(1);
        rhsIdentiferVector.reserveItems(size);
        passiveValueVector.reserveItems(size);

        indexManager.get().assignIndex(lhsIndex);
        Real& primalEntry = primals[lhsIndex];
        cast().pushStmtData(lhsIndex, size, primalEntry, manualStatementReversal);

        primalEntry = lhsValue;
      }

      /*******************************************************************************
       * Section: Function from PositionalEvaluationTapeInterface
       *
       */

      CODI_INLINE void evaluate(Position const& start, Position const& end) {
        checkAdjointSize(indexManager.get().getLargestAssignedIndex());

        evaluate(start, end, adjoints.data());
      }

      CODI_INLINE void resetTo(Position const& pos) {

        internalResetPrimalValues(pos);

        Base::resetTo(pos);
      }


      /*******************************************************************************
       * Section: Function from PreaccumulationEvaluationTapeInterface
       *
       */

      void evaluateKeepState(Position const& start, Position const& end) {
        checkAdjointSize(indexManager.get().getLargestAssignedIndex());

        internalEvaluateReverse<false>(start, end, adjoints.data());

        if(!TapeTypes::IsLinearIndexHandler) {

          internalEvaluatePrimal(end, start);
        }
      }

      void evaluateForwardKeepState(Position const& start, Position const& end) {
        checkAdjointSize(indexManager.get().getLargestAssignedIndex());

        if(!TapeTypes::IsLinearIndexHandler) {
          internalResetPrimalValues(end);
        }

        internalEvaluateForward<false>(start, end, adjoints.data());
      }

      /*******************************************************************************
       * Section: Function from PrimalEvaluationTapeInterface
       *
       */

    protected:

      template<typename Rhs>
      static Real statementPrimal(
          Real* primalVector, Config::ArgumentSize numberOfPassiveArguments,
          size_t& curConstantPos, PassiveReal const* const constantValues,
          size_t& curPassivePos, Real const* const passiveValues,
          size_t& curRhsIdentifiersPos, Identifier const* const rhsIdentifiers) {

        using Construtor = ConstructStaticContextlLogic<Rhs, Impl, 0, 0>;
        using StaticRhs = typename Construtor::ResultType;

        size_t constexpr MaxActiveArgs = MaxNumberOfActiveTypeArguments<Rhs>::value;
        size_t constexpr MaxConstantArgs = MaxNumberOfConstantArguments<Rhs>::value;

        for(Config::ArgumentSize curPos = 0; curPos < numberOfPassiveArguments; curPos += 1) {
          primalVector[curPos] = passiveValues[curPassivePos + curPos];
        }

        StaticRhs staticsRhs = Construtor::construct(
              primalVector,
              &rhsIdentifiers[curRhsIdentifiersPos],
              &constantValues[curConstantPos]);

        // Adapt vector positions
        curConstantPos += MaxConstantArgs;
        curPassivePos += numberOfPassiveArguments;
        curRhsIdentifiersPos += MaxActiveArgs;

        return staticsRhs.getValue();
      }

      WRAP_FUNCTION(Wrap_internalEvaluatePrimalStack, Impl::internalEvaluatePrimalStack);

      CODI_INLINE static void internalEvaluatePrimalVector(NestedPosition const& start, NestedPosition const& end,
                                                         Real* primalData,
                                                         ConstantValueVector& constantValueVector) {
        Wrap_internalEvaluatePrimalStack evalFunc{};
        constantValueVector.evaluateForward(start, end, evalFunc, primalData);
      }

      WRAP_FUNCTION(Wrap_internalEvaluatePrimalVector, internalEvaluatePrimalVector);

    public:

      using Base::evaluatePrimal;
      CODI_NO_INLINE void evaluatePrimal(Position const& start, Position const& end) {

        // TODO: implement primal value only accessor
        PrimalAdjointVectorAccess<Real, Identifier, Gradient> primalAdjointAccess(adjoints.data(), primals);

        if(TapeTypes::IsLinearIndexHandler) {

          Wrap_internalEvaluatePrimalVector evalFunc{};
          Base::internalEvaluateExtFuncPrimal(start, end, evalFunc,
                                              &primalAdjointAccess, primals.data(), constantValueVector);
        } else {
          Base::internalEvaluateExtFuncPrimal(start, end, PrimalValueBaseTape::internalEvaluatePrimalVector,
                                              &primalAdjointAccess, primals.data(), constantValueVector);
        }
      }

      Real& primal(Identifier const& identifier) {

        return primals[identifier];
      }

      Real const& primal(Identifier const& identifier) const {

        return primals[identifier];
      }

    private:

      CODI_INLINE void checkAdjointSize(Identifier const& identifier) {
        if(identifier >= (Identifier)adjoints.size()) {
          resizeAdjointsVector();
        }
      }

      CODI_INLINE void checkPrimalSize(bool generatedNewIndex) {
        if(TapeTypes::IsLinearIndexHandler) {
          if(indexManager.get().getLargestAssignedIndex() >= (Identifier)primals.size()) {
            resizePrimalVector(primals.size() + Config::ChunkSize);
          }
        } else {
          if(generatedNewIndex) {
            resizePrimalVector(indexManager.get().getLargestAssignedIndex() + 1);
          }
        }
      }

      CODI_NO_INLINE void resizeAdjointsVector() {
        adjoints.resize(indexManager.get().getLargestAssignedIndex() + 1);
      }

      CODI_NO_INLINE void resizePrimalVector(size_t newSize) {
        primals.resize(newSize);
      }
  };
}
