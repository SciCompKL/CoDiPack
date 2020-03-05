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
#include "data/chunk.hpp"
#include "data/chunkVector.hpp"
#include "indices/indexManagerInterface.hpp"
#include "interfaces/reverseTapeInterface.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  template<typename _Real, typename _Gradient, typename _IndexManager>
  struct PrimalValueTapeTypes {
    public:

      using Real = DECLARE_DEFAULT(_Real, double);
      using Gradient = DECLARE_DEFAULT(_Gradient, double);
      using IndexManager = DECLARE_DEFAULT(_IndexManager, TEMPLATE(IndexManagerInterface<int>));

      using Identifier = typename IndexManager::Index;
      using PassiveReal = PassiveRealType<Real>;

      constexpr static bool IsLinearIndexHandler = IndexManager::IsLinear;
      constexpr static bool IsStaticIndexHandler = !IsLinearIndexHandler;

      using EvalPointer = void (*)(
          Real* primalVector, Gradient* adjointVector,
          Gradient lhsAdjoint, Config::ArgumentSize numberOfPassiveArguments,
          size_t& curConstantPos, PassiveReal const* const constantValues,
          size_t& curPassivePos, Real const* const passiveValues,
          size_t& curRhsIdentifiersPos, Identifier const* const rhsIdentifiers);

      using StatementChunk = Chunk4<Identifier, Config::ArgumentSize, Real, EvalPointer>;
      using StatementVector = ChunkVector<StatementChunk, IndexManager>;

      using IdentifierChunk = Chunk1<Identifier>;
      using RhsIdentifierVector = ChunkVector<IdentifierChunk, StatementVector>;

      using PassiveValueChunk = Chunk1<Real>;
      using PassiveValueVector = ChunkVector<PassiveValueChunk, RhsIdentifierVector>;

      using ConstantValueChunk = Chunk1<PassiveReal>;
      using ConstantValueVector = ChunkVector<ConstantValueChunk, PassiveValueVector>;
  };

  template<typename _TapeTypes>
  struct PrimalValueReuseTape :
      public ReverseTapeInterface<
          typename _TapeTypes::Real,
          typename _TapeTypes::Gradient,
          typename _TapeTypes::Identifier>
  {
    public:

      using TapeTypes = DECLARE_DEFAULT(_TapeTypes, TEMPLATE(PrimalValueTapeTypes<double, double, IndexManagerInterface<int>));

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

      static bool constexpr AllowJacobianOptimization = false;

    protected:

      MemberStore<IndexManager, PrimalValueReuseTape, TapeTypes::IsStaticIndexHandler> indexManager;
      StatementVector statementVector;
      RhsIdentifierVector rhsIdentiferVector;
      PassiveValueVector passiveValueVector;
      ConstantValueVector constantValueVector;

      bool active;

      std::vector<Gradient> adjoints;
      std::vector<Real> primals;

    public:

      PrimalValueReuseTape() :
        indexManager(Config::MaxArgumentSize),
        statementVector(Config::ChunkSize),
        rhsIdentiferVector(Config::ChunkSize),
        passiveValueVector(Config::ChunkSize),
        constantValueVector(Config::ChunkSize),
        active(false),
        adjoints(1),
        primals()
      {
        resizePrimalVector();

        statementVector.setNested(&indexManager.get());
        rhsIdentiferVector.setNested(&statementVector);
        passiveValueVector.setNested(&rhsIdentiferVector);
        constantValueVector.setNested(&passiveValueVector);
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

      template<typename Lhs, typename Rhs>
      CODI_INLINE void store(LhsExpressionInterface<Real, Gradient, PrimalValueReuseTape, Lhs>& lhs,
                 ExpressionInterface<Real, Rhs> const& rhs) {

        if(active) {
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
            if(generatedNewIndex) {
              resizePrimalVector();
            }

            Real& primalEntry = primals[lhs.cast().getIdentifier()];
            statementVector.pushData(lhs.cast().getIdentifier(), passiveArguments, primalEntry, statementReversal<Rhs>);

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
      CODI_INLINE void store(LhsExpressionInterface<Real, Gradient, PrimalValueReuseTape, Lhs>& lhs,
                 LhsExpressionInterface<Real, Gradient, PrimalValueReuseTape, Rhs> const& rhs) {

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
      CODI_INLINE void store(LhsExpressionInterface<Real, Gradient, PrimalValueReuseTape, Lhs>& lhs, PassiveReal const& rhs) {
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

      struct IncrementReversalLogic : public TraversalLogic<IncrementReversalLogic> {
        public:
          template<typename Node>
          CODI_INLINE enableIfStaticContextActiveType<Node> term(Node const& node, Real jacobian, Gradient const& lhsAdjoint, Gradient* adjointVector) {
            using std::isfinite;
            ENABLE_CHECK(Config::CheckZeroIndex, (Identifier)Config::MaxArgumentSize <= node.getIdentifier()) {
              ENABLE_CHECK(Config::IgnoreInvalidJacobies, isfinite(jacobian)) {
                ENABLE_CHECK(Config::CheckJacobiIsZero, !isTotalZero(jacobian)) {
                  adjointVector[node.getIdentifier()] += jacobian * lhsAdjoint;
                }
              }
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
          Real* primalVector, Gradient* adjointVector,
          Gradient lhsAdjoint, Config::ArgumentSize numberOfPassiveArguments,
          size_t& curConstantPos, PassiveReal const* const constantValues,
          size_t& curPassivePos, Real const* const passiveValues,
          size_t& curRhsIdentifiersPos, Identifier const* const rhsIdentifiers) {

        using Construtor = ConstructStaticContextlLogic<Rhs, PrimalValueReuseTape, 0, 0>;
        using StaticRhs = typename Construtor::ResultType;

        size_t constexpr MaxActiveArgs = MaxNumberOfActiveTypeArguments<Rhs>::value;
        size_t constexpr MaxConstantArgs = MaxNumberOfConstantArguments<Rhs>::value;

        // Adapt vector positions
        curConstantPos -= MaxConstantArgs;
        curPassivePos -= numberOfPassiveArguments;
        curRhsIdentifiersPos -= MaxActiveArgs;

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

      static void internalEvaluateReverse(
          /* data from call */
          Real* primalVector, Gradient* adjointVector,
          /* data constant value vector */
          size_t& curConstantPos, size_t const& endConstantPos, PassiveReal const* const constantValues,
          /* data passive value vector */
          size_t& curPassivePos, size_t const& endPassivePos, Real const* const passiveValues,
          /* data rhs identifiers vector */
          size_t& curRhsIdentifiersPos, size_t const& endRhsIdentifiersPos, Identifier const* const rhsIdentifiers,
          /* data statement vector */
          size_t& curStatementPos, size_t const& endStatementPos,
              Identifier const* const lhsIdentifiers,
              Config::ArgumentSize const* const numberOfPassiveArguments,
              Real const * const oldPrimalValues,
              EvalPointer const * const stmtEvalFunc
          ) {

        CODI_UNUSED(endConstantPos, endPassivePos, endRhsIdentifiersPos);

        while(curStatementPos > endStatementPos) {
          curStatementPos -= 1;

          Identifier const lhsIdentifier = lhsIdentifiers[curStatementPos];

          Gradient const lhsAdjoint = adjointVector[lhsIdentifier];
          adjointVector[lhsIdentifier] = Gradient();

          primalVector[lhsIdentifier] = oldPrimalValues[curStatementPos];

          stmtEvalFunc[curStatementPos](primalVector, adjointVector, lhsAdjoint,
              numberOfPassiveArguments[curStatementPos], curConstantPos, constantValues,
              curPassivePos, passiveValues, curRhsIdentifiersPos, rhsIdentifiers);
        }
      }

      CODI_NO_INLINE void evaluate() {
        checkAdjointSize(indexManager.get().getLargestAssignedIndex());

        std::vector primalsCopy(primals);

        constantValueVector.evaluateReverse(constantValueVector.getPosition(),
                                            constantValueVector.getZeroPosition(),
                                            internalEvaluateReverse,
                                            primalsCopy.data(), adjoints.data());
      }

      template<typename Lhs>
      void registerInput(LhsExpressionInterface<Real, Gradient, PrimalValueReuseTape, Lhs>& value) {
        bool generatedNewIndex = indexManager.get().assignUnusedIndex(value.cast().getIdentifier());

        if(generatedNewIndex) {
          resizePrimalVector();
        }

        primals[value.cast().getIdentifier()] = value.cast().value();
      }

      template<typename Lhs>
      CODI_INLINE void registerOutput(LhsExpressionInterface<Real, Gradient, PrimalValueReuseTape, Lhs>& value) {
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

        for(Real& primal : primals) {
          primal = Real();
        }

        constantValueVector.reset();
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

    private:

      CODI_INLINE void checkAdjointSize(Identifier const& identifier) {
        if(identifier >= (Identifier)adjoints.size()) {
          resizeAdjointsVector();
        }
      }

      CODI_NO_INLINE void resizeAdjointsVector() {
        adjoints.resize(indexManager.get().getLargestAssignedIndex() + 1);
      }

      CODI_NO_INLINE void resizePrimalVector() {
        primals.resize(indexManager.get().getLargestAssignedIndex() + 1);
      }
  };
}
