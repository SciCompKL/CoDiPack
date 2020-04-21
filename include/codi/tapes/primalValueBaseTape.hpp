#pragma once

#include <algorithm>
#include <functional>
#include <type_traits>

#include "../aux/macros.h"
#include "../aux/memberStore.hpp"
#include "../config.h"
#include "../expressions/lhsExpressionInterface.hpp"
#include "../expressions/logic/compileTimeTraversalLogic.hpp"
#include "../expressions/logic/helpers/forEachTermLogic.hpp"
#include "../expressions/logic/helpers/jacobianComputationLogic.hpp"
#include "../expressions/logic/traversalLogic.hpp"
#include "../expressions/logic/constructStaticContext.hpp"
#include "../traits/expressionTraits.hpp"
#include "aux/primalAdjointVectorAccess.hpp"
#include "data/chunk.hpp"
#include "indices/indexManagerInterface.hpp"
#include "commonTapeImplementation.hpp"
#include "statementEvaluators/statementEvaluatorTapeInterface.hpp"
#include "statementEvaluators/statementEvaluatorInterface.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  template<typename _Real, typename _Gradient, typename _IndexManager, template <typename> class _StatementEvaluator, template<typename, typename> class _Vector>
  struct PrimalValueTapeTypes : public TapeTypesInterface {
    public:

      using Real = DECLARE_DEFAULT(_Real, double);
      using Gradient = DECLARE_DEFAULT(_Gradient, double);
      using IndexManager = DECLARE_DEFAULT(_IndexManager, TEMPLATE(IndexManagerInterface<int>));
      using StatementEvaluator = DECLARE_DEFAULT(TEMPLATE(_StatementEvaluator<Real>), TEMPLATE(StatementEvaluatorInterface<double>));
      template<typename Chunk, typename Nested>
      using Vector = DECLARE_DEFAULT(TEMPLATE(_Vector<Chunk, Nested>), TEMPLATE(DataInterface<Nested>));

      using Identifier = typename IndexManager::Index;
      using PassiveReal = PassiveRealType<Real>;

      constexpr static bool IsLinearIndexHandler = IndexManager::IsLinear;
      constexpr static bool IsStaticIndexHandler = !IsLinearIndexHandler;

      using EvalHandle = typename StatementEvaluator::Handle;

      using StatementChunk = typename std::conditional<
                                IsLinearIndexHandler,
                                Chunk2<Config::ArgumentSize, EvalHandle>,
                                Chunk4<Identifier, Config::ArgumentSize, Real, EvalHandle>
                              >::type;
      using StatementVector = Vector<StatementChunk, IndexManager>;

      using IdentifierChunk = Chunk1<Identifier>;
      using RhsIdentifierVector = Vector<IdentifierChunk, StatementVector>;

      using PassiveValueChunk = Chunk1<Real>;
      using PassiveValueVector = Vector<PassiveValueChunk, RhsIdentifierVector>;

      using ConstantValueChunk = Chunk1<PassiveReal>;
      using ConstantValueVector = Vector<ConstantValueChunk, PassiveValueVector>;

      using NestedVector = ConstantValueVector;
  };

  template<typename _TapeTypes, typename _Impl>
  struct PrimalValueBaseTape :
      public CommonTapeImplementation<_TapeTypes, _Impl>,
      public StatementEvaluatorTapeInterface<typename _TapeTypes::Real>,
      public StatementEvaluatorInnerTapeInterface<typename _TapeTypes::Real>
  {
    public:

      using TapeTypes = DECLARE_DEFAULT(_TapeTypes, TEMPLATE(PrimalValueTapeTypes<double, double, IndexManagerInterface<int>));
      using Impl = DECLARE_DEFAULT(_Impl, TEMPLATE(ReverseTapeInterface<double, double, int>));

      using Base = CommonTapeImplementation<TapeTypes, Impl>;
      friend Base;

      using Real = typename TapeTypes::Real;
      using Gradient = typename TapeTypes::Gradient;
      using IndexManager = typename TapeTypes::IndexManager;
      using StatementEvaluator = typename TapeTypes::StatementEvaluator;
      using Identifier = typename TapeTypes::Identifier;

      using EvalHandle = typename TapeTypes::EvalHandle;

      using StatementVector = typename TapeTypes::StatementVector;
      using RhsIdentifierVector = typename TapeTypes::RhsIdentifierVector;
      using PassiveValueVector = typename TapeTypes::PassiveValueVector;
      using ConstantValueVector = typename TapeTypes::ConstantValueVector;

      using PassiveReal = PassiveRealType<Real>;

      using NestedPosition = typename ConstantValueVector::Position;
      using Position = typename Base::Position;

      static bool constexpr AllowJacobianOptimization = false;
      static bool constexpr RequiresPrimalRestore = !TapeTypes::IsLinearIndexHandler;

    protected:

      static EvalHandle const jacobianExpressions[Config::MaxArgumentSize];

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
          EvalHandle evalHandle);

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
        Base::options.insert(ConfigurationOption::LargestIdentifier);
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

      struct CountActiveArguments : public ForEachTermLogic<CountActiveArguments> {
        public:
          template<typename Node>
          CODI_INLINE void handleActive(Node const& node, size_t& numberOfActiveArguments) {
            ENABLE_CHECK(Config::CheckZeroIndex, 0 != node.getIdentifier()) {

              numberOfActiveArguments += 1;
            }
          }
      };

      struct PushIdentfierPassiveAndConstant : public ForEachTermLogic<PushIdentfierPassiveAndConstant> {
        public:
          template<typename Node>
          CODI_INLINE void handleActive(
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
          CODI_INLINE void handleConstant(
              Node const& node,
              RhsIdentifierVector& rhsIdentiferVector,
              PassiveValueVector& passiveValueVector,
              ConstantValueVector& constantValueVector,
              size_t& curPassiveArgument) {

            CODI_UNUSED(rhsIdentiferVector, passiveValueVector, curPassiveArgument);

            constantValueVector.pushData(node.getValue());
          }
      };

    public:

      template<typename Lhs, typename Rhs>
      CODI_INLINE void store(LhsExpressionInterface<Real, Gradient, Impl, Lhs>& lhs,
                 ExpressionInterface<Real, Rhs> const& rhs) {

        ENABLE_CHECK(Config::CheckTapeActivity, cast().isActive()) {
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

      template<typename Lhs, typename Rhs>
      CODI_INLINE void store(LhsExpressionInterface<Real, Gradient, Impl, Lhs>& lhs,
                 LhsExpressionInterface<Real, Gradient, Impl, Rhs> const& rhs) {

        ENABLE_CHECK(Config::CheckTapeActivity, cast().isActive()) {
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
          cast().pushStmtData(value.cast().getIdentifier(), Config::StatementInputTag, primalEntry,
                              StatementEvaluator::template createHandle<Impl, Impl, Lhs>());
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



      struct IncrementReversalLogic : public JacobianComputationLogic<Real, IncrementReversalLogic> {
        public:
          template<typename Node>
          CODI_INLINE void handleJacobianOnActive(Node const& node, Real jacobian, Gradient const& lhsAdjoint, ADJOINT_VECTOR_TYPE* adjointVector) {
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
      };

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

      struct IncrementForwardLogic : public JacobianComputationLogic<Real, IncrementForwardLogic> {
        public:
          template<typename Node>
          CODI_INLINE void handleJacobianOnActive(Node const& node, Real jacobian, Gradient& lhsTangent, ADJOINT_VECTOR_TYPE* adjointVector) {
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
      };

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
          case ConfigurationOption::LargestIdentifier:
            return indexManager.get().getLargestAssignedIndex();
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
            adjoints.resize(value);
            break;
          case ConfigurationOption::ConstantValuesSize:
            constantValueVector.resize(value);
            break;
          case ConfigurationOption::LargestIdentifier:
            CODI_EXCEPTION("Tried to set a get only option.");
            break;
          case ConfigurationOption::PassiveValuesSize:
            passiveValueVector.resize(value);
            break;
          case ConfigurationOption::RhsIdentifiersSize:
            rhsIdentiferVector.resize(value);
            break;
          case ConfigurationOption::PrimalSize:
            primals.resize(value);
            break;
          case ConfigurationOption::StatementSize:
            return statementVector.resize(value);
          default:
            Base::setOption(option, value);
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

      template<size_t _size>
      struct JacobianStatementGenerator {
        public:

          static size_t constexpr size = _size;

          template<typename Expr, typename ... Args>
          static Real statementEvaluateForward(Args&& ... args) {
            CODI_EXCEPTION("Forward evaluation of jacobian statement not possible.");
          }

          template<typename Expr, typename ... Args>
          static Real statementEvaluatePrimal(Args&& ... args) {
            CODI_EXCEPTION("Primal evaluation of jacobian statement not possible.");

          }

          template<typename Expr>
          static void statementEvaluateReverse(
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

          template<typename Expr, typename ... Args>
          static Real statementEvaluateForwardInner(Args&& ... args) {
            CODI_EXCEPTION("Forward evaluation of jacobian statement not possible.");
          }

          template<typename Expr, typename ... Args>
          static Real statementEvaluatePrimalInner(Args&& ... args) {
            CODI_EXCEPTION("Primal evaluation of jacobian statement not possible.");
          }

          template<typename Expr>
          static void statementEvaluateReverseInner(
              Real* primalVector, ADJOINT_VECTOR_TYPE* adjointVector,
              Gradient lhsAdjoint,
              size_t& curConstantPos, PassiveReal const* const constantValues,
              size_t& curRhsIdentifiersPos, Identifier const* const rhsIdentifiers) {

            CODI_UNUSED(primalVector, curConstantPos, constantValues);

            size_t curPrimalPos = size;
            size_t endPos = curRhsIdentifiersPos - size;

            #if CODI_VariableAdjointInterfaceInPrimalTapes
              bool const lhsZero = adjointVector->isLhsZero();
            #else
              bool const lhsZero = isTotalZero(lhsAdjoint);
            #endif

            ENABLE_CHECK(Config::SkipZeroAdjointEvaluation, !lhsZero) {
              while(curRhsIdentifiersPos > endPos) {
                curRhsIdentifiersPos -= 1;
                curPrimalPos -= 1;

                Real const& jacobian = primalVector[curPrimalPos];
                #if CODI_VariableAdjointInterfaceInPrimalTapes
                  adjointVector->updateAdjointWithLhs(rhsIdentifiers[curRhsIdentifiersPos], jacobian);
                #else
                  adjointVector[rhsIdentifiers[curRhsIdentifiersPos]] += jacobian * lhsAdjoint;
                #endif

              }
            } else {
              curRhsIdentifiersPos -= size;
            }

          }
      };

    public:

      void pushJacobiManual(Real const& jacobi, Real const& value, Identifier const& index) {
        CODI_UNUSED(value);

        passiveValueVector.pushData(jacobi);
        rhsIdentiferVector.pushData(index);
      }

      void storeManual(Real const& lhsValue, Identifier& lhsIndex, Config::ArgumentSize const& size) {
        CODI_UNUSED(lhsValue);

        statementVector.reserveItems(1);
        rhsIdentiferVector.reserveItems(size);
        passiveValueVector.reserveItems(size);

        indexManager.get().assignIndex(lhsIndex);
        Real& primalEntry = primals[lhsIndex];
        cast().pushStmtData(lhsIndex, size, primalEntry, PrimalValueBaseTape::jacobianExpressions[size]);

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

      /*******************************************************************************
       * Section: Function from StatementEvaluatorTapeInterface
       *
       */

      template<typename Rhs>
      static Real statementEvaluateForwardInner(
          Real* primalVector, ADJOINT_VECTOR_TYPE* adjointVector,
          Gradient& lhsTangent,
          size_t& curConstantPos, PassiveReal const* const constantValues,
          size_t& curRhsIdentifiersPos, Identifier const* const rhsIdentifiers) {

        using Construtor = ConstructStaticContextlLogic<Rhs, Impl, 0, 0>;
        using StaticRhs = typename Construtor::ResultType;

        StaticRhs staticsRhs = Construtor::construct(
              primalVector,
              &rhsIdentifiers[curRhsIdentifiersPos],
              &constantValues[curConstantPos]);

        IncrementForwardLogic incrementForward;

        incrementForward.eval(staticsRhs, 1.0, lhsTangent, adjointVector);
        return staticsRhs.getValue();

      }

      template<typename Func>
      static Real statementEvaluateForwardFull(
          Func const& evalInner, size_t const& maxActiveArgs, size_t const& maxConstantArgs,
          Real* primalVector, ADJOINT_VECTOR_TYPE* adjointVector,
          Gradient& lhsTangent, Config::ArgumentSize numberOfPassiveArguments,
          size_t& curConstantPos, PassiveReal const* const constantValues,
          size_t& curPassivePos, Real const* const passiveValues,
          size_t& curRhsIdentifiersPos, Identifier const* const rhsIdentifiers) {

        for(Config::ArgumentSize curPos = 0; curPos < numberOfPassiveArguments; curPos += 1) {
          primalVector[curPos] = passiveValues[curPassivePos + curPos];

        }

        Real ret = evalInner(primalVector, adjointVector, lhsTangent,
                             curConstantPos, constantValues,
                             curRhsIdentifiersPos, rhsIdentifiers);

        // Adapt vector positions
        curConstantPos += maxConstantArgs;
        curPassivePos += numberOfPassiveArguments;
        curRhsIdentifiersPos += maxActiveArgs;

        return ret;
      }

      template<typename Rhs>
      static Real statementEvaluateForward(
          Real* primalVector, ADJOINT_VECTOR_TYPE* adjointVector,
          Gradient& lhsTangent, Config::ArgumentSize numberOfPassiveArguments,
          size_t& curConstantPos, PassiveReal const* const constantValues,
          size_t& curPassivePos, Real const* const passiveValues,
          size_t& curRhsIdentifiersPos, Identifier const* const rhsIdentifiers) {

        size_t constexpr MaxActiveArgs = MaxNumberOfActiveTypeArguments<Rhs>::value;
        size_t constexpr MaxConstantArgs = MaxNumberOfConstantArguments<Rhs>::value;

        return statementEvaluateForwardFull(statementEvaluateForwardInner<Rhs>, MaxActiveArgs, MaxConstantArgs,
                                            primalVector, adjointVector, lhsTangent, numberOfPassiveArguments,
                                            curConstantPos, constantValues,
                                            curPassivePos, passiveValues,
                                            curRhsIdentifiersPos, rhsIdentifiers);
      }

      template<typename Rhs>
      static Real statementEvaluatePrimalInner(
          Real* primalVector,
          size_t& curConstantPos, PassiveReal const* const constantValues,
          size_t& curRhsIdentifiersPos, Identifier const* const rhsIdentifiers) {

        using Construtor = ConstructStaticContextlLogic<Rhs, Impl, 0, 0>;
        using StaticRhs = typename Construtor::ResultType;

        StaticRhs staticsRhs = Construtor::construct(
              primalVector,
              &rhsIdentifiers[curRhsIdentifiersPos],
              &constantValues[curConstantPos]);

        return staticsRhs.getValue();
      }

      template<typename Func>
      static Real statementEvaluatePrimalFull(
          Func const& evalInner, size_t const& maxActiveArgs, size_t const& maxConstantArgs,
          Real* primalVector, Config::ArgumentSize numberOfPassiveArguments,
          size_t& curConstantPos, PassiveReal const* const constantValues,
          size_t& curPassivePos, Real const* const passiveValues,
          size_t& curRhsIdentifiersPos, Identifier const* const rhsIdentifiers) {

        for(Config::ArgumentSize curPos = 0; curPos < numberOfPassiveArguments; curPos += 1) {
          primalVector[curPos] = passiveValues[curPassivePos + curPos];
        }

        Real ret = evalInner(primalVector, curConstantPos, constantValues, curRhsIdentifiersPos, rhsIdentifiers);

        // Adapt vector positions
        curConstantPos += maxConstantArgs;
        curPassivePos += numberOfPassiveArguments;
        curRhsIdentifiersPos += maxActiveArgs;

        return ret;
      }

      template<typename Rhs>
      static Real statementEvaluatePrimal(
          Real* primalVector, Config::ArgumentSize numberOfPassiveArguments,
          size_t& curConstantPos, PassiveReal const* const constantValues,
          size_t& curPassivePos, Real const* const passiveValues,
          size_t& curRhsIdentifiersPos, Identifier const* const rhsIdentifiers) {

        size_t constexpr MaxActiveArgs = MaxNumberOfActiveTypeArguments<Rhs>::value;
        size_t constexpr MaxConstantArgs = MaxNumberOfConstantArguments<Rhs>::value;

        return statementEvaluatePrimalFull(statementEvaluatePrimalInner<Rhs>, MaxActiveArgs, MaxConstantArgs,
                                    primalVector, numberOfPassiveArguments,
                                    curConstantPos, constantValues,
                                    curPassivePos, passiveValues,
                                    curRhsIdentifiersPos, rhsIdentifiers);
      }

      template<typename Rhs>
      CODI_INLINE static void statementEvaluateReverseInner(
          Real* primalVector, ADJOINT_VECTOR_TYPE* adjointVector,
          Gradient lhsAdjoint,
          size_t& curConstantPos, PassiveReal const* const constantValues,
          size_t& curRhsIdentifiersPos, Identifier const* const rhsIdentifiers) {

        using Construtor = ConstructStaticContextlLogic<Rhs, Impl, 0, 0>;
        using StaticRhs = typename Construtor::ResultType;

        StaticRhs staticsRhs = Construtor::construct(
              primalVector,
              &rhsIdentifiers[curRhsIdentifiersPos],
              &constantValues[curConstantPos]);

        IncrementReversalLogic incrementReverse;

        incrementReverse.eval(staticsRhs, 1.0, const_cast<Gradient const&>(lhsAdjoint), adjointVector);
      }


      template<typename Func>
      CODI_INLINE static void statementEvaluateReverseFull(
          Func const& evalInner, size_t const& maxActiveArgs, size_t const& maxConstantArgs,
          Real* primalVector, ADJOINT_VECTOR_TYPE* adjointVector,
          Gradient lhsAdjoint, Config::ArgumentSize numberOfPassiveArguments,
          size_t& curConstantPos, PassiveReal const* const constantValues,
          size_t& curPassivePos, Real const* const passiveValues,
          size_t& curRhsIdentifiersPos, Identifier const* const rhsIdentifiers) {

        // Adapt vector positions
        curConstantPos -= maxConstantArgs;
        curPassivePos -= numberOfPassiveArguments;
        curRhsIdentifiersPos -= maxActiveArgs;

        ENABLE_CHECK(Config::SkipZeroAdjointEvaluation, !isTotalZero(lhsAdjoint)) {
          for(Config::ArgumentSize curPos = 0; curPos < numberOfPassiveArguments; curPos += 1) {
            primalVector[curPos] = passiveValues[curPassivePos + curPos];
          }

          evalInner(primalVector, adjointVector, lhsAdjoint, curConstantPos, constantValues,
                                        curRhsIdentifiersPos, rhsIdentifiers);
        }
      }

      template<typename Rhs>
      CODI_INLINE static void statementEvaluateReverse(
          Real* primalVector, ADJOINT_VECTOR_TYPE* adjointVector,
          Gradient lhsAdjoint, Config::ArgumentSize numberOfPassiveArguments,
          size_t& curConstantPos, PassiveReal const* const constantValues,
          size_t& curPassivePos, Real const* const passiveValues,
          size_t& curRhsIdentifiersPos, Identifier const* const rhsIdentifiers) {

        size_t constexpr maxActiveArgs = MaxNumberOfActiveTypeArguments<Rhs>::value;
        size_t constexpr maxConstantArgs = MaxNumberOfConstantArguments<Rhs>::value;
        statementEvaluateReverseFull(
              statementEvaluateReverseInner<Rhs>, maxActiveArgs, maxConstantArgs,
              primalVector, adjointVector, lhsAdjoint, numberOfPassiveArguments,
              curConstantPos, constantValues,
              curPassivePos, passiveValues,
              curRhsIdentifiersPos, rhsIdentifiers);
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

  #define CREATE_EXPRESSION(size) \
      TapeTypes::StatementEvaluator::template createHandle<Impl, typename Impl::template JacobianStatementGenerator<size>, void>()

  template <typename TapeTypes, typename Impl>
  const typename TapeTypes::EvalHandle PrimalValueBaseTape<TapeTypes, Impl>::jacobianExpressions[Config::MaxArgumentSize] = {
    CREATE_EXPRESSION(0), CREATE_EXPRESSION(1), CREATE_EXPRESSION(2), CREATE_EXPRESSION(3), CREATE_EXPRESSION(4),
    CREATE_EXPRESSION(5), CREATE_EXPRESSION(6), CREATE_EXPRESSION(7), CREATE_EXPRESSION(8), CREATE_EXPRESSION(9),
    CREATE_EXPRESSION(10), CREATE_EXPRESSION(11), CREATE_EXPRESSION(12), CREATE_EXPRESSION(13), CREATE_EXPRESSION(14),
    CREATE_EXPRESSION(15), CREATE_EXPRESSION(16), CREATE_EXPRESSION(17), CREATE_EXPRESSION(18), CREATE_EXPRESSION(19),
    CREATE_EXPRESSION(20), CREATE_EXPRESSION(21), CREATE_EXPRESSION(22), CREATE_EXPRESSION(23), CREATE_EXPRESSION(24),
    CREATE_EXPRESSION(25), CREATE_EXPRESSION(26), CREATE_EXPRESSION(27), CREATE_EXPRESSION(28), CREATE_EXPRESSION(29),
    CREATE_EXPRESSION(30), CREATE_EXPRESSION(31), CREATE_EXPRESSION(32), CREATE_EXPRESSION(33), CREATE_EXPRESSION(34),
    CREATE_EXPRESSION(35), CREATE_EXPRESSION(36), CREATE_EXPRESSION(37), CREATE_EXPRESSION(38), CREATE_EXPRESSION(39),
    CREATE_EXPRESSION(40), CREATE_EXPRESSION(41), CREATE_EXPRESSION(42), CREATE_EXPRESSION(43), CREATE_EXPRESSION(44),
    CREATE_EXPRESSION(45), CREATE_EXPRESSION(46), CREATE_EXPRESSION(47), CREATE_EXPRESSION(48), CREATE_EXPRESSION(49),
    CREATE_EXPRESSION(50), CREATE_EXPRESSION(51), CREATE_EXPRESSION(52), CREATE_EXPRESSION(53), CREATE_EXPRESSION(54),
    CREATE_EXPRESSION(55), CREATE_EXPRESSION(56), CREATE_EXPRESSION(57), CREATE_EXPRESSION(58), CREATE_EXPRESSION(59),
    CREATE_EXPRESSION(60), CREATE_EXPRESSION(61), CREATE_EXPRESSION(62), CREATE_EXPRESSION(63), CREATE_EXPRESSION(64),
    CREATE_EXPRESSION(65), CREATE_EXPRESSION(66), CREATE_EXPRESSION(67), CREATE_EXPRESSION(68), CREATE_EXPRESSION(69),
    CREATE_EXPRESSION(70), CREATE_EXPRESSION(71), CREATE_EXPRESSION(72), CREATE_EXPRESSION(73), CREATE_EXPRESSION(74),
    CREATE_EXPRESSION(75), CREATE_EXPRESSION(76), CREATE_EXPRESSION(77), CREATE_EXPRESSION(78), CREATE_EXPRESSION(79),
    CREATE_EXPRESSION(80), CREATE_EXPRESSION(81), CREATE_EXPRESSION(82), CREATE_EXPRESSION(83), CREATE_EXPRESSION(84),
    CREATE_EXPRESSION(85), CREATE_EXPRESSION(86), CREATE_EXPRESSION(87), CREATE_EXPRESSION(88), CREATE_EXPRESSION(89),
    CREATE_EXPRESSION(90), CREATE_EXPRESSION(91), CREATE_EXPRESSION(92), CREATE_EXPRESSION(93), CREATE_EXPRESSION(94),
    CREATE_EXPRESSION(95), CREATE_EXPRESSION(96), CREATE_EXPRESSION(97), CREATE_EXPRESSION(98), CREATE_EXPRESSION(99),
    CREATE_EXPRESSION(100), CREATE_EXPRESSION(101), CREATE_EXPRESSION(102), CREATE_EXPRESSION(103), CREATE_EXPRESSION(104),
    CREATE_EXPRESSION(105), CREATE_EXPRESSION(106), CREATE_EXPRESSION(107), CREATE_EXPRESSION(108), CREATE_EXPRESSION(109),
    CREATE_EXPRESSION(110), CREATE_EXPRESSION(111), CREATE_EXPRESSION(112), CREATE_EXPRESSION(113), CREATE_EXPRESSION(114),
    CREATE_EXPRESSION(115), CREATE_EXPRESSION(116), CREATE_EXPRESSION(117), CREATE_EXPRESSION(118), CREATE_EXPRESSION(119),
    CREATE_EXPRESSION(120), CREATE_EXPRESSION(121), CREATE_EXPRESSION(122), CREATE_EXPRESSION(123), CREATE_EXPRESSION(124),
    CREATE_EXPRESSION(125), CREATE_EXPRESSION(126), CREATE_EXPRESSION(127), CREATE_EXPRESSION(128), CREATE_EXPRESSION(129),
    CREATE_EXPRESSION(130), CREATE_EXPRESSION(131), CREATE_EXPRESSION(132), CREATE_EXPRESSION(133), CREATE_EXPRESSION(134),
    CREATE_EXPRESSION(135), CREATE_EXPRESSION(136), CREATE_EXPRESSION(137), CREATE_EXPRESSION(138), CREATE_EXPRESSION(139),
    CREATE_EXPRESSION(140), CREATE_EXPRESSION(141), CREATE_EXPRESSION(142), CREATE_EXPRESSION(143), CREATE_EXPRESSION(144),
    CREATE_EXPRESSION(145), CREATE_EXPRESSION(146), CREATE_EXPRESSION(147), CREATE_EXPRESSION(148), CREATE_EXPRESSION(149),
    CREATE_EXPRESSION(150), CREATE_EXPRESSION(151), CREATE_EXPRESSION(152), CREATE_EXPRESSION(153), CREATE_EXPRESSION(154),
    CREATE_EXPRESSION(155), CREATE_EXPRESSION(156), CREATE_EXPRESSION(157), CREATE_EXPRESSION(158), CREATE_EXPRESSION(159),
    CREATE_EXPRESSION(160), CREATE_EXPRESSION(161), CREATE_EXPRESSION(162), CREATE_EXPRESSION(163), CREATE_EXPRESSION(164),
    CREATE_EXPRESSION(165), CREATE_EXPRESSION(166), CREATE_EXPRESSION(167), CREATE_EXPRESSION(168), CREATE_EXPRESSION(169),
    CREATE_EXPRESSION(170), CREATE_EXPRESSION(171), CREATE_EXPRESSION(172), CREATE_EXPRESSION(173), CREATE_EXPRESSION(174),
    CREATE_EXPRESSION(175), CREATE_EXPRESSION(176), CREATE_EXPRESSION(177), CREATE_EXPRESSION(178), CREATE_EXPRESSION(179),
    CREATE_EXPRESSION(180), CREATE_EXPRESSION(181), CREATE_EXPRESSION(182), CREATE_EXPRESSION(183), CREATE_EXPRESSION(184),
    CREATE_EXPRESSION(185), CREATE_EXPRESSION(186), CREATE_EXPRESSION(187), CREATE_EXPRESSION(188), CREATE_EXPRESSION(189),
    CREATE_EXPRESSION(190), CREATE_EXPRESSION(191), CREATE_EXPRESSION(192), CREATE_EXPRESSION(193), CREATE_EXPRESSION(194),
    CREATE_EXPRESSION(195), CREATE_EXPRESSION(196), CREATE_EXPRESSION(197), CREATE_EXPRESSION(198), CREATE_EXPRESSION(199),
    CREATE_EXPRESSION(200), CREATE_EXPRESSION(201), CREATE_EXPRESSION(202), CREATE_EXPRESSION(203), CREATE_EXPRESSION(204),
    CREATE_EXPRESSION(205), CREATE_EXPRESSION(206), CREATE_EXPRESSION(207), CREATE_EXPRESSION(208), CREATE_EXPRESSION(209),
    CREATE_EXPRESSION(210), CREATE_EXPRESSION(211), CREATE_EXPRESSION(212), CREATE_EXPRESSION(213), CREATE_EXPRESSION(214),
    CREATE_EXPRESSION(215), CREATE_EXPRESSION(216), CREATE_EXPRESSION(217), CREATE_EXPRESSION(218), CREATE_EXPRESSION(219),
    CREATE_EXPRESSION(220), CREATE_EXPRESSION(221), CREATE_EXPRESSION(222), CREATE_EXPRESSION(223), CREATE_EXPRESSION(224),
    CREATE_EXPRESSION(225), CREATE_EXPRESSION(226), CREATE_EXPRESSION(227), CREATE_EXPRESSION(228), CREATE_EXPRESSION(229),
    CREATE_EXPRESSION(230), CREATE_EXPRESSION(231), CREATE_EXPRESSION(232), CREATE_EXPRESSION(233), CREATE_EXPRESSION(234),
    CREATE_EXPRESSION(235), CREATE_EXPRESSION(236), CREATE_EXPRESSION(237), CREATE_EXPRESSION(238), CREATE_EXPRESSION(239),
    CREATE_EXPRESSION(240), CREATE_EXPRESSION(241), CREATE_EXPRESSION(242), CREATE_EXPRESSION(243), CREATE_EXPRESSION(244),
    CREATE_EXPRESSION(245), CREATE_EXPRESSION(246), CREATE_EXPRESSION(247), CREATE_EXPRESSION(248), CREATE_EXPRESSION(249),
    CREATE_EXPRESSION(250), CREATE_EXPRESSION(251), CREATE_EXPRESSION(252), CREATE_EXPRESSION(253)
  };

  #undef CREATE_EXPRESSION
}
