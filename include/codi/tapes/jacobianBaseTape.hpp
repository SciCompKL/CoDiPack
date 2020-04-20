#pragma once

#include <algorithm>
#include <type_traits>

#include "../aux/macros.h"
#include "../aux/memberStore.hpp"
#include "../config.h"
#include "../expressions/lhsExpressionInterface.hpp"
#include "../expressions/logic/compileTimeTraversalLogic.hpp"
#include "../expressions/logic/traversalLogic.hpp"
#include "../expressions/logic/helpers/forEachTermLogic.hpp"
#include "../expressions/logic/helpers/jacobianComputationLogic.hpp"
#include "../expressions/referenceActiveType.hpp"
#include "../traits/expressionTraits.hpp"
#include "aux/adjointVectorAccess.hpp"
#include "aux/jacobianSorter.hpp"
#include "data/chunk.hpp"
#include "indices/indexManagerInterface.hpp"
#include "commonTapeImplementation.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  template<typename _Real, typename _Gradient, typename _IndexManager, template<typename, typename> class _Vector>
  struct JacobianTapeTypes : public TapeTypesInterface {
    public:

      using Real = DECLARE_DEFAULT(_Real, double);
      using Gradient = DECLARE_DEFAULT(_Gradient, double);
      using IndexManager = DECLARE_DEFAULT(_IndexManager, TEMPLATE(IndexManagerInterface<int>));
      template<typename Chunk, typename Nested>
      using Vector = DECLARE_DEFAULT(TEMPLATE(_Vector<Chunk, Nested>), TEMPLATE(DataInterface<Nested>));

      using Identifier = typename IndexManager::Index;

      constexpr static bool IsLinearIndexHandler = IndexManager::IsLinear;
      constexpr static bool IsStaticIndexHandler = !IsLinearIndexHandler;

      using StatementChunk = typename std::conditional<
                                 IsLinearIndexHandler,
                                 Chunk1<Config::ArgumentSize>,
                                 Chunk2<Identifier, Config::ArgumentSize>
                               >::type;
      using StatementVector = Vector<StatementChunk, IndexManager>;

      using JacobianChunk = Chunk2<Real, Identifier>;
      using JacobianVector = Vector<JacobianChunk, StatementVector>;

      using NestedVector = JacobianVector;
  };

  template<typename _TapeTypes, typename _Impl>
  struct JacobianBaseTape : public CommonTapeImplementation<_TapeTypes, _Impl> {
    public:

      using TapeTypes = DECLARE_DEFAULT(_TapeTypes, TEMPLATE(JacobianTapeTypes<double, double, IndexManagerInterface<int>));
      using Impl = DECLARE_DEFAULT(_Impl, TEMPLATE(FullTapeInterface<double, double, int, EmptyPosition>));

      using Base = CommonTapeImplementation<TapeTypes, Impl>;
      friend Base;

      using Real = typename TapeTypes::Real;
      using Gradient = typename TapeTypes::Gradient;
      using IndexManager = typename TapeTypes::IndexManager;
      using Identifier = typename TapeTypes::Identifier;

      using StatementVector = typename TapeTypes::StatementVector;
      using JacobianVector = typename TapeTypes::JacobianVector;

      using PassiveReal = PassiveRealType<Real>;

      using NestedPosition = typename JacobianVector::Position;
      using Position = typename Base::Position;

      static bool constexpr AllowJacobianOptimization = true;

    protected:


#if CODI_CombineJacobianArguments
      JacobianSorter<Real, Identifier> jacobianSorter;
#endif

      MemberStore<IndexManager, Impl, TapeTypes::IsStaticIndexHandler> indexManager;
      StatementVector statementVector;
      JacobianVector jacobianVector;

      std::vector<Gradient> adjoints;

    public:

      CODI_INLINE Impl const& cast() const {
        return static_cast<Impl const&>(*this);
      }

      CODI_INLINE Impl& cast() {
        return static_cast<Impl&>(*this);
      }

    protected:

      /*******************************************************************************
       * Section: Methods expected in the child class.
       *
       * Description: TODO
       *
       */

      template<typename ... Args>
      static void internalEvaluateForward(Args&& ... args);

      template<typename ... Args>
      static void internalEvaluateReverse(Args&& ... args);

      void pushStmtData(Identifier const& index, Config::ArgumentSize const& numberOfArguments);

    public:

      JacobianBaseTape() :
        Base(),
  #if CODI_CombineJacobianArguments
        jacobianSorter(),
  #endif
        indexManager(0),
        statementVector(Config::ChunkSize),
        jacobianVector(Config::ChunkSize),
        adjoints(1)
      {
        statementVector.setNested(&indexManager.get());
        jacobianVector.setNested(&statementVector);

        Base::init(&jacobianVector);

        Base::options.insert(ConfigurationOption::AdjointSize);
        Base::options.insert(ConfigurationOption::JacobianSize);
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

      struct PushJacobianLogic : public JacobianComputationLogic<Real, PushJacobianLogic> {
        public:
          template<typename Node, typename DataVector>
          CODI_INLINE void handleJacobianOnActive(Node const& node, Real jacobian, DataVector& dataVector) {
            using std::isfinite;
            ENABLE_CHECK(Config::CheckZeroIndex, 0 != node.getIdentifier()) {
              ENABLE_CHECK(Config::IgnoreInvalidJacobies, isfinite(jacobian)) {
                ENABLE_CHECK(Config::CheckJacobiIsZero, !isTotalZero(jacobian)) {
                  dataVector.pushData(jacobian, node.getIdentifier());
                }
              }
            }
          }

          template<typename Type, typename DataVector>
          CODI_INLINE void handleJacobianOnActive(ReferenceActiveType<Type> const& node, Real jacobian, DataVector& dataVector) {
            CODI_UNUSED(dataVector);

            using std::isfinite;
            ENABLE_CHECK(Config::IgnoreInvalidJacobies, isfinite(jacobian)) {
              // Do a delayed push for these termination nodes, accumulate the jacobian in the local member
              node.jacobian += jacobian;
            }
          }
      };

      struct PushDelayedJacobianLogic : public ForEachTermLogic<PushDelayedJacobianLogic> {
        public:
          template<typename Type, typename DataVector>
          CODI_INLINE void handleActive(ReferenceActiveType<Type> const& node, DataVector& dataVector) {
            using std::isfinite;
            ENABLE_CHECK(Config::CheckZeroIndex, 0 != node.getIdentifier()) {
              ENABLE_CHECK(Config::CheckJacobiIsZero, !isTotalZero(node.jacobian)) {
                dataVector.pushData(node.jacobian, node.getIdentifier());

                // Reset the jacobian here such that it is not pushed multiple times and ready for the next store
                node.jacobian = Real();
              }
            }
          }

          using ForEachTermLogic<PushDelayedJacobianLogic>::handleActive;
      };

      template<typename Rhs>
      CODI_INLINE void pushJacobians(ExpressionInterface<Real, Rhs> const& rhs) {
        PushJacobianLogic pushJacobianLogic;
        PushDelayedJacobianLogic pushDelayedJacobianLogic;

#if CODI_CombineJacobianArguments
        auto& insertVector = jacobianSorter;
#else
        auto& insertVector = jacobianVector;
#endif

        pushJacobianLogic.eval(rhs.cast(), 1.0, insertVector);
        pushDelayedJacobianLogic.eval(rhs.cast(), insertVector);

#if CODI_CombineJacobianArguments
        jacobianSorter.storeData(jacobianVector);
#endif
      }

    public:

      template<typename Lhs, typename Rhs>
      CODI_INLINE void store(LhsExpressionInterface<Real, Gradient, Impl, Lhs>& lhs,
                 ExpressionInterface<Real, Rhs> const& rhs) {

        ENABLE_CHECK(Config::CheckTapeActivity, cast().isActive()) {
          size_t constexpr MaxArgs = MaxNumberOfActiveTypeArguments<Rhs>::value;

          statementVector.reserveItems(1);
          typename JacobianVector::InternalPosHandle jacobianStart = jacobianVector.reserveItems(MaxArgs);

          pushJacobians(rhs);

          size_t numberOfArguments = jacobianVector.getPushedDataCount(jacobianStart);
          if(0 != numberOfArguments) {
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

      template<typename Lhs>
      CODI_INLINE void registerInput(LhsExpressionInterface<Real, Gradient, Impl, Lhs>& value) {
        indexManager.get().assignUnusedIndex(value.cast().getIdentifier());

        if(TapeTypes::IsLinearIndexHandler) {
          statementVector.reserveItems(1);
          cast().pushStmtData(value.cast().getIdentifier(), Config::StatementInputTag);
        }
      }

      CODI_INLINE void clearAdjoints() {
        for(Gradient& gradient : adjoints) {
          gradient = Gradient();
        }
      }

    protected:

      CODI_INLINE TapeValues internalGetTapeValues() const {
        std::string name;
        if(TapeTypes::IsLinearIndexHandler) {
          name = "CoDi Tape Statistics ( JacobiLinearTape )";
        } else {
          name = "CoDi Tape Statistics ( JacobiReuseTape )";
        }
        TapeValues values = TapeValues(name);

        size_t nAdjoints      = indexManager.get().getLargestAssignedIndex();
        double memoryAdjoints = static_cast<double>(nAdjoints) * static_cast<double>(sizeof(Gradient)) * TapeValues::BYTE_TO_MB;

        values.addSection("Adjoint vector");
        values.addUnsignedLongEntry("Number of adjoints", nAdjoints);
        values.addDoubleEntry("Memory allocated", memoryAdjoints, true, true);

        values.addSection("Index manager");
        indexManager.get().addToTapeValues(values);

        values.addSection("Statement entries");
        statementVector.addToTapeValues(values);
        values.addSection("Jacobian entries");
        jacobianVector.addToTapeValues(values);

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
      CODI_INLINE static void incrementAdjoints(
          Adjoint* adjointVector,
          Adjoint const& lhsAdjoint,
          Config::ArgumentSize const& numberOfArguments,
          size_t& curJacobianPos,
          Real const* const rhsJacobians,
          Identifier const* const rhsIdentifiers) {

        size_t endJacobianPos = curJacobianPos - numberOfArguments;

        ENABLE_CHECK(Config::SkipZeroAdjointEvaluation, !isTotalZero(lhsAdjoint)){
          while(endJacobianPos < curJacobianPos) {
            curJacobianPos -= 1;
            adjointVector[rhsIdentifiers[curJacobianPos]] += rhsJacobians[curJacobianPos] * lhsAdjoint;
          }
        } else {
          curJacobianPos = endJacobianPos;
        }
      }

    public:
      template<typename Adjoint>
      CODI_NO_INLINE void evaluate(Position const& start, Position const& end, Adjoint* data) {
        AdjointVectorAccess<Real, Identifier, Adjoint> adjointWrapper(data);

        auto evalFunc = [this] (NestedPosition const& start, NestedPosition const& end,
                                Adjoint* data) {
            jacobianVector.evaluateReverse(start, end, Impl::template internalEvaluateReverse<Adjoint>, data);
        };
        Base::internalEvaluateExtFunc(start, end, evalFunc, &adjointWrapper, data);
      }

    protected:

      template<typename Adjoint>
      CODI_INLINE static void incrementTangents (
          Adjoint const* const adjointVector,
          Adjoint& lhsAdjoint,
          Config::ArgumentSize const& numberOfArguments,
          size_t& curJacobianPos,
          Real const* const rhsJacobians,
          Identifier const* const rhsIdentifiers) {

        size_t endJacobianPos = curJacobianPos + numberOfArguments;

        while(curJacobianPos < endJacobianPos) {
          curJacobianPos += 1;
          lhsAdjoint += rhsJacobians[curJacobianPos] * adjointVector[rhsIdentifiers[curJacobianPos]];
        }
      }

    public:
      template<typename Adjoint>
      CODI_NO_INLINE void evaluateForward(Position const& start, Position const& end, Adjoint* data) {
        AdjointVectorAccess<Real, Identifier, Adjoint> adjointWrapper(data);

        auto evalFunc = [this] (NestedPosition const& start, NestedPosition const& end, Adjoint* data) {
          jacobianVector.evaluateForward(start, end, Impl::template internalEvaluateForward<Adjoint>, data);
        };
        Base::internalEvaluateExtFunc(start, end, evalFunc, &adjointWrapper, data, jacobianVector);

      }

      /*******************************************************************************
       * Section: Function from DataManagementTapeInterface
       *
       */

      CODI_INLINE void swap(Impl& other) {

        // Index manager does not need to be swapped, it is either static or swapped in with the vector data
        // Vectors are swapped recursively in the base class

        std::swap(adjoints, other.adjoints);

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
          case ConfigurationOption::JacobianSize:
            return jacobianVector.getDataSize();
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
          case ConfigurationOption::JacobianSize:
            return jacobianVector.resize(value);
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
        cast().registerInput(value);

        return Real();
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

      void pushJacobiManual(Real const& jacobi, Real const& value, Gradient const& index) {
        CODI_UNUSED(value);

        jacobianVector.pushData(jacobi, index);
      }

      void storeManual(Real const& lhsValue, Gradient& lhsIndex, Config::ArgumentSize const& size) {
        CODI_UNUSED(lhsValue);

        statementVector.reserveItems(1);
        jacobianVector.reserveItems(size);

        indexManager.get().assignIndex(lhsIndex);
        cast().pushStmtData(lhsIndex, (Config::ArgumentSize)size);
      }

    public:

      /*******************************************************************************
       * Section: Function from PositionalEvaluationTapeInterface
       *
       */

      CODI_INLINE void evaluate(Position const& start, Position const& end) {
        checkAdjointSize(indexManager.get().getLargestAssignedIndex());

        evaluate(start, end, adjoints.data());
      }

      /*******************************************************************************
       * Section: Function from PreaccumulationEvaluationTapeInterface
       *
       */

      void evaluateKeepState(Position const& start, Position const& end) {
        evaluate(start, end);
      }

      void evaluateForwardKeepState(Position const& start, Position const& end) {
        evaluateForward(start, end);
      }

      /*******************************************************************************
       * Section: Function from PrimalEvaluationTapeInterface
       *
       */

      using Base::evaluatePrimal;
      void evaluatePrimal(Position const& start, Position const& end) {
        CODI_UNUSED(start, end);

        CODI_EXCEPTION("Accessing primal evaluation of an Jacobian tape.");
      }

      Real& primal(Identifier const& identifier) {
        CODI_UNUSED(identifier);

        CODI_EXCEPTION("Accessing primal vector of an Jacobian tape.");

        static Real temp;
        return temp;
      }

      Real const& primal(Identifier const& identifier) const {
        CODI_UNUSED(identifier);

        CODI_EXCEPTION("Accessing primal vector of an Jacobian tape.");

        return Real();
      }

    private:

      CODI_INLINE void checkAdjointSize(Identifier const& identifier) {
        if(identifier >= (Identifier)adjoints.size()) {
          resizeAdjointsVector();
        }
      }

      CODI_NO_INLINE void resizeAdjointsVector() {
        adjoints.resize(indexManager.get().getLargestAssignedIndex() + 1);
      }
  };
}

