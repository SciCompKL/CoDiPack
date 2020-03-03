#pragma once

#include <algorithm>
#include <type_traits>

#include "../aux/macros.h"
#include "../aux/memberStore.hpp"
#include "../config.h"
#include "../expressions/lhsExpressionInterface.hpp"
#include "../expressions/logic/compileTimeTraversalLogic.hpp"
#include "../expressions/logic/traversalLogic.hpp"
#include "../traits/expressionTraits.hpp"
#include "data/chunk.hpp"
#include "data/chunkVector.hpp"
#include "indices/indexManagerInterface.hpp"
#include "interfaces/reverseTapeInterface.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  template<typename _Real, typename _Gradient, typename _IndexManager>
  struct JacobianTapeTypes {
    public:

      using Real = DECLARE_DEFAULT(_Real, double);
      using Gradient = DECLARE_DEFAULT(_Gradient, double);
      using IndexManager = DECLARE_DEFAULT(_IndexManager, TEMPLATE(IndexManagerInterface<int>));

      using Identifier = typename IndexManager::Index;

      constexpr static bool IsLinearIndexHandler = IndexManager::IsLinear;
      constexpr static bool IsStaticIndexHandler = !IsLinearIndexHandler;

      using StatementChunk = typename std::conditional<
                                 IsLinearIndexHandler,
                                 Chunk1<Config::ArgumentSize>,
                                 Chunk2<Identifier, Config::ArgumentSize>
                               >::type;
      using StatementVector = ChunkVector<StatementChunk, IndexManager>;

      using JacobianChunk = Chunk2<Real, Identifier>;
      using JacobianVector = ChunkVector<JacobianChunk, StatementVector>;

  };

  template<typename _TapeTypes, typename _Impl>
  struct JacobianBaseTape :
      public ReverseTapeInterface<
          typename _TapeTypes::Real,
          typename _TapeTypes::Gradient,
          typename _TapeTypes::Identifier>
  {
    public:

      using TapeTypes = DECLARE_DEFAULT(_TapeTypes, TEMPLATE(JacobianTapeTypes<double, double, IndexManagerInterface<int>));
      using Impl = DECLARE_DEFAULT(_Impl, TEMPLATE(ReverseTapeInterface<double, double, int>));

      using Real = typename TapeTypes::Real;
      using Gradient = typename TapeTypes::Gradient;
      using IndexManager = typename TapeTypes::IndexManager;
      using Identifier = typename TapeTypes::Identifier;

      using StatementVector = typename TapeTypes::StatementVector;
      using JacobianVector = typename TapeTypes::JacobianVector;

      using PassiveReal = PassiveRealType<Real>;

      static bool constexpr AllowJacobianOptimization = true;

    protected:

      MemberStore<IndexManager, Impl, TapeTypes::IsStaticIndexHandler> indexManager;
      StatementVector statementVector;
      JacobianVector jacobianVector;

      bool active;

      std::vector<Gradient> adjoints;

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

    public:
      template<typename Lhs>
      void registerInput(LhsExpressionInterface<Real, Gradient, Impl, Lhs>& value);

    protected:

      void pushStmtData(Identifier const& index, Config::ArgumentSize const& numberOfArguments);

      template<typename ... Args>
      static void internalEvaluateRevere(Args&& ... args);

    public:

      JacobianBaseTape() :
        indexManager(0),
        statementVector(Config::ChunkSize),
        jacobianVector(Config::ChunkSize),
        active(false),
        adjoints(1)
      {
        statementVector.setNested(&indexManager.get());
        jacobianVector.setNested(&statementVector);
      }

      CODI_INLINE void setGradient(Identifier const& identifier, Gradient const& gradient) {
        gradient(identifier) = gradient;
      }

      CODI_INLINE Gradient const& getGradient(Identifier const& identifier) const {
        return gradient(identifier);
      }

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

      struct PushJacobianLogic : public TraversalLogic<PushJacobianLogic> {
        public:
          template<typename Node>
          CODI_INLINE enableIfLhsExpression<Real, Gradient, Impl, Node> term(Node const& node, Real jacobian, JacobianVector& jacobianVector, size_t& numberOfArguments) {
            using std::isfinite;
            ENABLE_CHECK(Config::CheckZeroIndex, 0 != node.getIdentifier()) {
              ENABLE_CHECK(Config::IgnoreInvalidJacobies, isfinite(jacobian)) {
                ENABLE_CHECK(Config::CheckJacobiIsZero, !isTotalZero(jacobian)) {
                  jacobianVector.pushData(jacobian, node.getIdentifier());

                  numberOfArguments += 1;
                }
              }
            }
          }
          using TraversalLogic<PushJacobianLogic>::term;

          template<size_t LeafNumber, typename Leaf, typename Root>
          CODI_INLINE void link(Leaf const& leaf, Root const& root, Real jacobian, JacobianVector& jacobianVector, size_t& numberOfArguments) {

            Real curJacobian = root.template getJacobian<LeafNumber>() * jacobian;

            this->toNode(leaf, curJacobian, jacobianVector, numberOfArguments);
          }
      };

      struct MaxNumberOfArguments : public CompileTimeTraversalLogic<size_t, MaxNumberOfArguments> {
        public:
          template<typename Node, typename = enableIfLhsExpression<Real, Gradient, Impl, Node>>
          CODI_INLINE static constexpr size_t term() {
            return 1;
          }
          using CompileTimeTraversalLogic<size_t, MaxNumberOfArguments>::term;
      };

      template<typename Lhs, typename Rhs>
      CODI_INLINE void store(LhsExpressionInterface<Real, Gradient, Impl, Lhs>& lhs,
                 ExpressionInterface<Real, Rhs> const& rhs) {

        if(active) {
          PushJacobianLogic pushJacobianLogic;
          size_t constexpr MaxArgs = MaxNumberOfArguments::template eval<Rhs>();

          statementVector.reserveItems(1);
          jacobianVector.reserveItems(MaxArgs);

          size_t numberOfArguments = 0;
          pushJacobianLogic.eval(rhs.cast(), 1.0, jacobianVector, numberOfArguments);

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

        if(active) {
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

      CODI_INLINE static void incrementAdjoints(
          Gradient* adjointVector,
          Gradient const& lhsAdjoint,
          Config::ArgumentSize const& numberOfArguments,
          size_t& curJacobianPos,
          Real const* const rhsJacobians,
          Identifier const* const rhsIdentifiers) {


        ENABLE_CHECK(Config::SkipZeroAdjointEvaluation, !isTotalZero(lhsAdjoint)){
          for(Config::ArgumentSize argPos = 0; argPos < numberOfArguments; argPos += 1) {
            curJacobianPos -= 1;
            adjointVector[rhsIdentifiers[curJacobianPos]] += rhsJacobians[curJacobianPos] * lhsAdjoint;
          }
        } else {
          curJacobianPos -= numberOfArguments;
        }
      }

      CODI_NO_INLINE void evaluate() {
        checkAdjointSize(indexManager.get().getLargestAssignedIndex());

        jacobianVector.evaluateReverse(jacobianVector.getPosition(),
                                       jacobianVector.getZeroPosition(),
                                       Impl::internalEvaluateRevere,
                                       adjoints.data());
      }

      template<typename Lhs>
      CODI_INLINE void registerOutput(LhsExpressionInterface<Real, Gradient, Impl, Lhs>& value) {
        store<Lhs, Lhs>(value, static_cast<ExpressionInterface<Real, Lhs> const&>(value));
      }

      CODI_INLINE void setActive() {
        active = true;
      }

      CODI_INLINE void setPassive() {
        active = false;
      }

      CODI_INLINE bool isActive() const {
        return active;
      }

      void clearAdjoints() {
        for(Gradient& gradient : adjoints) {
          gradient = Gradient();
        }
      }

      void reset(bool resetAdjoints = true) {
        if(resetAdjoints) {
          clearAdjoints();
        }

        jacobianVector.reset();

        // Requires extra reset since the default vector implementation forwards to resetTo
        indexManager.get().reset();
      }

      template<typename Stream = std::ostream> void printStatistics(Stream& out = std::cout) const {
        getTapeValues().formatDefault(out);
      }

      template<typename Stream = std::ostream> void printTableHeader(Stream& out = std::cout) const {
        getTapeValues().formatHeader(out);
      }

      template<typename Stream = std::ostream> void printTableRow(Stream& out = std::cout) const {
        getTapeValues().formatRow(out);
      }

      TapeValues getTapeValues() const {
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

