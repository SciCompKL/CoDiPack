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
#include "indices/indexManagerInterface.hpp"
#include "primalValueBaseTape.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  template<typename _TapeTypes>
  struct PrimalValueReuseTape : public PrimalValueBaseTape<_TapeTypes, PrimalValueReuseTape<_TapeTypes>> {
    public:

      using TapeTypes = DECLARE_DEFAULT(_TapeTypes, TEMPLATE(PrimalValueTapeTypes<double, double, IndexManagerInterface<int>));
      using Base = PrimalValueBaseTape<_TapeTypes, PrimalValueReuseTape<_TapeTypes>>;
      friend Base;

      using Real = typename TapeTypes::Real;
      using Gradient = typename TapeTypes::Gradient;
      using Identifier = typename TapeTypes::Identifier;
      using PassiveReal = PassiveRealType<Real>;
      using StatementEvaluator = typename TapeTypes::StatementEvaluator;
      using EvalHandle = typename TapeTypes::EvalHandle;
      using Position = typename Base::Position;

      using StatementVector = typename TapeTypes::StatementVector;

      PrimalValueReuseTape() : Base() {}

      using Base::clearAdjoints;
      void clearAdjoints(Position const& start, Position const& end) {

        // clear adjoints
        auto clearFunc = [this] (Identifier* lhsIndex, Config::ArgumentSize* passiveArgs, Real* oldPrimal, EvalHandle* evalHandle) {
          CODI_UNUSED(passiveArgs, oldPrimal, evalHandle);

          if(*lhsIndex < this->adjoints.size()) {
            this->adjoints[*lhsIndex] = Gradient();
          }
        };

        using StmtPosition = typename StatementVector::Position;
        StmtPosition startStmt = this->externalFunctionVector.template extractPosition<StmtPosition>(start);
        StmtPosition endStmt = this->externalFunctionVector.template extractPosition<StmtPosition>(end);

        this->statementVector.forEachReverse(startStmt, endStmt, clearFunc);
      }

    protected:

      CODI_INLINE static void internalEvaluateForwardStack(
          /* data from call */
          Real* primalVector, ADJOINT_VECTOR_TYPE* adjointVector,
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
              Real * const oldPrimalValues,
              EvalHandle const * const stmtEvalhandle
          ) {

        CODI_UNUSED(endConstantPos, endPassivePos, endRhsIdentifiersPos);

        while(curStatementPos < endStatementPos) {

          Identifier const lhsIdentifier = lhsIdentifiers[curStatementPos];


          Gradient lhsTangent = Gradient();

          oldPrimalValues[curStatementPos] = primalVector[lhsIdentifier];
          primalVector[lhsIdentifier] = StatementEvaluator::template callForward<PrimalValueReuseTape>(
                stmtEvalhandle[curStatementPos], primalVector, adjointVector, lhsTangent,
                numberOfPassiveArguments[curStatementPos], curConstantPos, constantValues,
                curPassivePos, passiveValues, curRhsIdentifiersPos, rhsIdentifiers);

          #if CODI_VariableAdjointInterfaceInPrimalTapes
            adjointVector->setLhsTangent(lhsIdentifier);
          #else
            adjointVector[lhsIdentifier] = lhsTangent;
          #endif
        }

        curStatementPos += 1;
      }

      CODI_INLINE static void internalEvaluatePrimalStack(
          /* data from call */
          Real* primalVector,
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
              Real * const oldPrimalValues,
              EvalHandle const * const stmtEvalhandle
          ) {

        CODI_UNUSED(endConstantPos, endPassivePos, endRhsIdentifiersPos);

        while(curStatementPos < endStatementPos) {

          Identifier const lhsIdentifier = lhsIdentifiers[curStatementPos];


          oldPrimalValues[curStatementPos] = primalVector[lhsIdentifier];
          primalVector[lhsIdentifier] = StatementEvaluator::template callPrimal<PrimalValueReuseTape>(
                stmtEvalhandle[curStatementPos], primalVector,
                numberOfPassiveArguments[curStatementPos], curConstantPos, constantValues,
                curPassivePos, passiveValues, curRhsIdentifiersPos, rhsIdentifiers);
        }

        curStatementPos += 1;
      }

      CODI_INLINE static void internalEvaluateReverseStack(
          /* data from call */
          Real* primalVector, ADJOINT_VECTOR_TYPE* adjointVector,
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
              EvalHandle const * const stmtEvalhandle
          ) {

        CODI_UNUSED(endConstantPos, endPassivePos, endRhsIdentifiersPos);

        while(curStatementPos > endStatementPos) {
          curStatementPos -= 1;

          Identifier const lhsIdentifier = lhsIdentifiers[curStatementPos];


          #if CODI_VariableAdjointInterfaceInPrimalTapes
            Gradient const lhsAdjoint;
            adjointVector->setLhsAdjoint(lhsIdentifier);
          #else
            Gradient const lhsAdjoint = adjointVector[lhsIdentifier];
            adjointVector[lhsIdentifier] = Gradient();
          #endif

          primalVector[lhsIdentifier] = oldPrimalValues[curStatementPos];

          StatementEvaluator::template callReverse<PrimalValueReuseTape>(
                stmtEvalhandle[curStatementPos], primalVector, adjointVector, lhsAdjoint,
                numberOfPassiveArguments[curStatementPos], curConstantPos, constantValues,
                curPassivePos, passiveValues, curRhsIdentifiersPos, rhsIdentifiers);
        }
      }

      CODI_INLINE void internalResetPrimalValues(Position const& pos) {

        // reset primals
        auto clearFunc = [this] (Identifier* lhsIndex, Config::ArgumentSize* passiveArgs, Real* oldPrimal, EvalHandle* evalHandle) {
          CODI_UNUSED(passiveArgs, evalHandle);

          this->primals[*lhsIndex] = *oldPrimal;
        };

        using StmtPosition = typename StatementVector::Position;
        StmtPosition startStmt = this->externalFunctionVector.template extractPosition<StmtPosition>(this->getPosition());
        StmtPosition endStmt = this->externalFunctionVector.template extractPosition<StmtPosition>(pos);

        this->statementVector.forEachReverse(startStmt, endStmt, clearFunc);

      }

      void pushStmtData(
          Identifier const& index,
          Config::ArgumentSize const& numberOfPassiveArguments,
          Real const& oldPrimalValue,
          EvalHandle evalHandle)
      {
        Base::statementVector.pushData(index, numberOfPassiveArguments, oldPrimalValue, evalHandle);
      }
  };
}
