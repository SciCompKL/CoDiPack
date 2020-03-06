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
  struct PrimalValueReuseTape : public PrimalValueBaseTape<_TapeTypes, PrimalValueReuseTape<_TapeTypes>> {
    public:

      using TapeTypes = DECLARE_DEFAULT(_TapeTypes, TEMPLATE(PrimalValueTapeTypes<double, double, IndexManagerInterface<int>));
      using Base = PrimalValueBaseTape<_TapeTypes, PrimalValueReuseTape<_TapeTypes>>;
      friend Base;

      using Real = typename TapeTypes::Real;
      using Gradient = typename TapeTypes::Gradient;
      using Identifier = typename TapeTypes::Identifier;
      using PassiveReal = PassiveRealType<Real>;
      using EvalPointer = typename TapeTypes::EvalPointer;

      PrimalValueReuseTape() : Base() {}

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

      void pushStmtData(
          Identifier const& index,
          Config::ArgumentSize const& numberOfPassiveArguments,
          Real const& oldPrimalValue,
          EvalPointer evalFunc)
      {
        Base::statementVector.pushData(index, numberOfPassiveArguments, oldPrimalValue, evalFunc);
      }
  };
}
