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
      using Identifier = typename TapeTypes::Identifier;
      using PassiveReal = PassiveRealType<Real>;
      using EvalPointer = typename TapeTypes::EvalPointer;

      PrimalValueLinearTape() : Base() {}

    protected:

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

          Gradient const lhsAdjoint = adjointVector[curAdjointPos];

          if(Config::StatementInputTag != nPassiveValues) {
            adjointVector[curAdjointPos] = Gradient();

            stmtEvalFunc[curStatementPos](primalVector, adjointVector, lhsAdjoint,
                nPassiveValues, curConstantPos, constantValues,
                curPassivePos, passiveValues, curRhsIdentifiersPos, rhsIdentifiers);
          }

          curAdjointPos -= 1;
        }
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
