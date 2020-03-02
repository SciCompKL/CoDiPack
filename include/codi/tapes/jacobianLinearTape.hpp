#pragma once

#include <algorithm>
#include <type_traits>

#include "../aux/macros.h"
#include "../config.h"
#include "../expressions/lhsExpressionInterface.hpp"
#include "../expressions/logic/compileTimeTraversalLogic.hpp"
#include "../expressions/logic/traversalLogic.hpp"
#include "../traits/expressionTraits.hpp"
#include "data/chunk.hpp"
#include "data/chunkVector.hpp"
#include "indices/linearIndexManager.hpp"
#include "interfaces/reverseTapeInterface.hpp"

#include "jacobianTape.hpp"

/** \copydoc codi::Namespace */
namespace codi {


  template<typename _Real, typename _Gradient, typename _Index>
  struct JacobianTape<_Real, _Gradient, LinearIndexManager<_Index>> : ReverseTapeInterface<_Real, _Gradient, _Index> {
    public:

      using Real = DECLARE_DEFAULT(_Real, double);
      using Gradient = DECLARE_DEFAULT(_Gradient, double);
      using IndexManager = DECLARE_DEFAULT(LinearIndexManager<_Index>, TEMPLATE(IndexManagerInterface<int>));

      using Identifier = typename IndexManager::Index;
      using PassiveReal = PassiveRealType<Real>;

      using StatementChunk = Chunk1<Config::ArgumentSize>;
      using StatementVector = ChunkVector<StatementChunk, IndexManager>;

      using JacobianChunk = Chunk2<Real, Identifier>;
      using JacobianVector = ChunkVector<JacobianChunk, StatementVector>;

      static bool constexpr AllowJacobianOptimization = true;

    private:

      IndexManager indexManager;

      StatementVector statementVector;
      JacobianVector jacobianVector;

      bool active;

      std::vector<Gradient> adjoints;

    public:

      JacobianTape() :
        indexManager(0),
        statementVector(Config::SmallChunkSize),
        jacobianVector(Config::ChunkSize),
        active(false),
        adjoints(1)
      {
        statementVector.setNested(&indexManager);
        jacobianVector.setNested(&statementVector);
      }

      void setGradient(Identifier const& identifier, Gradient const& gradient) {
        gradient(identifier) = gradient;
      }

      Gradient const& getGradient(Identifier const& identifier) const {
        return gradient(identifier);
      }

      Gradient& gradient(Identifier const& identifier) {
        checkAdjointSize(identifier);

        return adjoints[identifier];
      }
      Gradient const& gradient(Identifier const& identifier) const {
        if(identifier > (Identifier)adjoints.size()) {
          return adjoints[0];
        } else {
          return adjoints[identifier];
        }
      }

      template<typename Real>
      void initIdentifier(Real& value, Identifier& identifier) {
        CODI_UNUSED(value);

        identifier = IndexManager::UnusedIndex;
      }

      template<typename Real>
      void destroyIdentifier(Real& value, Identifier& identifier) {
        CODI_UNUSED(value);

        indexManager.freeIndex(identifier);
      }

      struct PushJacobianLogic : public TraversalLogic<PushJacobianLogic> {
        public:
          template<typename Node>
          CODI_INLINE enableIfLhsExpression<Real, Gradient, JacobianTape, Node> term(Node const& node, Real jacobian, JacobianVector& jacobianVector, size_t& numberOfArguments) {
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
          template<typename Node, typename = enableIfLhsExpression<Real, Gradient, JacobianTape, Node>>
          CODI_INLINE static constexpr size_t term() {
            return 1;
          }
          using CompileTimeTraversalLogic<size_t, MaxNumberOfArguments>::term;
      };

      template<typename Lhs, typename Rhs>
      void store(LhsExpressionInterface<Real, Gradient, JacobianTape, Lhs>& lhs,
                 ExpressionInterface<Real, Rhs> const& rhs) {

        if(active) {
          PushJacobianLogic pushJacobianLogic;
          size_t constexpr MaxArgs = MaxNumberOfArguments::template eval<Rhs>();

          statementVector.reserveItems(1);
          jacobianVector.reserveItems(MaxArgs);

          size_t numberOfArguments = 0;
          pushJacobianLogic.eval(rhs.cast(), 1.0, jacobianVector, numberOfArguments);

          if(0 != numberOfArguments) {
            indexManager.assignIndex(lhs.cast().getIdentifier());
            statementVector.pushData((Config::ArgumentSize)numberOfArguments);
          } else {
            indexManager.freeIndex(lhs.cast().getIdentifier());
          }
        } else {
          indexManager.freeIndex(lhs.cast().getIdentifier());
        }

        lhs.cast().value() = rhs.cast().getValue();
      }

      template<typename Lhs, typename Rhs>
      void store(LhsExpressionInterface<Real, Gradient, JacobianTape, Lhs>& lhs,
                 LhsExpressionInterface<Real, Gradient, JacobianTape, Rhs> const& rhs) {

        if(active) {
          if(IndexManager::AssignNeedsStatement || !Config::AssignOptimization) {
            store<Lhs, Rhs>(lhs, static_cast<ExpressionInterface<Real, Rhs> const&>(rhs));
          } else {
            indexManager.copyIndex(lhs.cast().getIdentifier(), rhs.cast().getIdentifier());
          }
        } else {
          indexManager.freeIndex(lhs.cast().getIdentifier());
        }

        lhs.cast().value() = rhs.cast().getValue();
      }

      template<typename Lhs>
      void store(LhsExpressionInterface<Real, Gradient, JacobianTape, Lhs>& lhs, PassiveReal const& rhs) {
        indexManager.freeIndex(lhs.cast().getIdentifier());

        lhs.cast().value() = rhs;
      }

      void evaluate() {
        checkAdjointSize(indexManager.getLargestAssignedIndex());

        auto evalReverse = [](
              /* data from call */
              Gradient* adjointVector,
              /* data from jacobian vector */
              size_t& curJacobianPos, size_t const& endJacobianPos, Real const* const rhsJacobians, Identifier const* const rhsIdentifiers ,
              /* data from statement vector */
              size_t& curStmtPos, size_t const& endStmtPos, Config::ArgumentSize const* const numberOfJacobians,
              /* data from index handler */
              size_t const& startAdjointPos, size_t const& endAdjointPos
            ) {

          CODI_UNUSED(endJacobianPos, endStmtPos);

          size_t curAdjointPos = startAdjointPos;

          while(curAdjointPos > endAdjointPos) {

            curStmtPos -= 1;
            Config::ArgumentSize const argsSize = numberOfJacobians[curStmtPos];

            Gradient const lhsAdjoint = adjointVector[curAdjointPos]; // Adjoint positions are shifted since we do not use the zero index

            if(Config::StatementInputTag != argsSize) {
              // No input value, perform regular statement evaluation
              adjointVector[curAdjointPos] = Gradient();

              curJacobianPos -= argsSize;

              ENABLE_CHECK(Config::SkipZeroAdjointEvaluation, !isTotalZero(lhsAdjoint)){
                for(Config::ArgumentSize argPos = 0; argPos < argsSize; argPos += 1) {
                  size_t curOffset = curJacobianPos + argPos;
                  adjointVector[rhsIdentifiers[curOffset]] += rhsJacobians[curOffset] * lhsAdjoint;
                }
              }
            }

            curAdjointPos -= 1;
          }
        };

        jacobianVector.evaluateReverse(jacobianVector.getPosition(),
                                       jacobianVector.getZeroPosition(),
                                       evalReverse,
                                       adjoints.data());
      }

      template<typename Lhs> void registerInput(LhsExpressionInterface<Real, Gradient, JacobianTape, Lhs>& value) {
        indexManager.assignUnusedIndex(value.cast().getIdentifier());
        statementVector.reserveItems(1);
        statementVector.pushData(Config::StatementInputTag);
      }

      template<typename Lhs> void registerOutput(LhsExpressionInterface<Real, Gradient, JacobianTape, Lhs>& value) {
        store<Lhs, Lhs>(value, static_cast<ExpressionInterface<Real, Lhs> const&>(value));
      }

      void setActive() {
        active = true;
      }

      void setPassive() {
        active = false;
      }

      bool isActive() const {
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
      }

      template<typename Stream = std::ostream> void printStatistics(Stream& out = std::cout) const { /*TODO: Implement */}
      template<typename Stream = std::ostream> void printTableHeader(Stream& out = std::cout) const { /*TODO: Implement */}
      template<typename Stream = std::ostream> void printTableRow(Stream& out = std::cout) const { /*TODO: Implement */}
      TapeValues getTapeValues() const { return TapeValues(); /*TODO: Implement */}

    private:

      CODI_INLINE void checkAdjointSize(Identifier const& identifier) {
        if(identifier >= (Identifier)adjoints.size()) {
          resizeAdjointsVector();
        }
      }

      CODI_NO_INLINE void resizeAdjointsVector() {
        adjoints.resize(indexManager.getLargestAssignedIndex() + 1);
      }
  };
}

