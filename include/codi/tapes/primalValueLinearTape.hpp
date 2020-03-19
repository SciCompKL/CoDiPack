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
#include "primalValueBaseTape.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  template<typename _TapeTypes>
  struct PrimalValueLinearTape : public PrimalValueBaseTape<_TapeTypes, PrimalValueLinearTape<_TapeTypes>> {
    public:

      using TapeTypes = DECLARE_DEFAULT(_TapeTypes, TEMPLATE(PrimalValueTapeTypes<double, double, IndexManagerInterface<int>));
      using Base = PrimalValueBaseTape<_TapeTypes, PrimalValueLinearTape<_TapeTypes>>;
      friend Base;

      using Real = typename TapeTypes::Real;
      using Gradient = typename TapeTypes::Gradient;
      using IndexManager = typename TapeTypes::IndexManager;
      using Identifier = typename TapeTypes::Identifier;
      using PassiveReal = PassiveRealType<Real>;
      using EvalPointer = typename TapeTypes::EvalPointer;
      using Position = typename Base::Position;

      PrimalValueLinearTape() : Base() {}

      using Base::clearAdjoints;
      void clearAdjoints(Position const& start, Position const& end) {

        using IndexPosition = typename IndexManager::Position;
        IndexPosition startIndex = this->externalFunctionVector.template extractPosition<IndexPosition>(start);
        IndexPosition endIndex = this->externalFunctionVector.template extractPosition<IndexPosition>(end);

        startIndex = std::min(startIndex, (IndexPosition)this->adjoints.size() - 1);
        endIndex = std::min(endIndex, (IndexPosition)this->adjoints.size() - 1);

        for(IndexPosition curPos = endIndex + 1; curPos <= startIndex; curPos += 1) {
          this->adjoints[curPos] = Gradient();
        }
      }

    protected:

      static void internalEvaluateForward(
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
              Config::ArgumentSize const* const numberOfPassiveArguments,
              EvalPointer const * const stmtEvalFunc,
          /* data from index handler */
          size_t const& startAdjointPos, size_t const& endAdjointPos
          ) {

        CODI_UNUSED(endConstantPos, endPassivePos, endRhsIdentifiersPos, endStatementPos);

        size_t curAdjointPos = startAdjointPos;

        while(curAdjointPos < endAdjointPos) {

          curAdjointPos += 1;

          Config::ArgumentSize nPassiveValues = numberOfPassiveArguments[curStatementPos];

          if(Config::StatementInputTag != nPassiveValues) {
            Gradient lhsTangent = Gradient();

            stmtEvalFunc[curStatementPos](primalVector, adjointVector, lhsTangent,
                nPassiveValues, curConstantPos, constantValues,
                curPassivePos, passiveValues, curRhsIdentifiersPos, rhsIdentifiers);

            #if CODI_VariableAdjointInterfaceInPrimalTapes
              adjointVector->setLhsTangent(curAdjointPos);
            #else
              adjointVector[curAdjointPos] = lhsTangent;
            #endif
          }

          curStatementPos += 1;
        }
      }

      static void internalEvaluateReverse(
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
              Config::ArgumentSize const* const numberOfPassiveArguments,
              EvalPointer const * const stmtEvalFunc,
          /* data from index handler */
          size_t const& startAdjointPos, size_t const& endAdjointPos
          ) {

        CODI_UNUSED(endConstantPos, endPassivePos, endRhsIdentifiersPos, endStatementPos);

        size_t curAdjointPos = startAdjointPos;

        while(curAdjointPos > endAdjointPos) {
          curStatementPos -= 1;

          Config::ArgumentSize nPassiveValues = numberOfPassiveArguments[curStatementPos];

          if(Config::StatementInputTag != nPassiveValues) {
            #if CODI_VariableAdjointInterfaceInPrimalTapes
              Gradient const lhsAdjoint;
              adjointVector->setLhsAdjoint(curAdjointPos);
            #else
              Gradient const lhsAdjoint = adjointVector[curAdjointPos];
              adjointVector[curAdjointPos] = Gradient();
            #endif

            stmtEvalFunc[curStatementPos](primalVector, adjointVector, lhsAdjoint,
                nPassiveValues, curConstantPos, constantValues,
                curPassivePos, passiveValues, curRhsIdentifiersPos, rhsIdentifiers);
          }

          curAdjointPos -= 1;
        }
      }

      void internalResetPrimalValues(Position const& pos) {
        CODI_UNUSED(pos);

        // Nothing to do
      }

      void pushStmtData(
          Identifier const& index,
          Config::ArgumentSize const& numberOfPassiveArguments,
          Real const& oldPrimalValue,
          EvalPointer evalFunc)
      {
        CODI_UNUSED(index, oldPrimalValue);

        Base::statementVector.pushData(numberOfPassiveArguments, evalFunc);
      }
  };
}
