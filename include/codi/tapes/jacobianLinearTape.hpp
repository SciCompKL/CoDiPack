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
#include "jacobianBaseTape.hpp"


/** \copydoc codi::Namespace */
namespace codi {


  template<typename _TapeTypes>
  struct JacobianLinearTape : public JacobianBaseTape<_TapeTypes, JacobianLinearTape<_TapeTypes>> {
    public:

      using TapeTypes = DECLARE_DEFAULT(_TapeTypes, TEMPLATE(JacobianTapeTypes<double, double, IndexManagerInterface<int>));

      using Base = JacobianBaseTape<_TapeTypes, JacobianLinearTape>;
      friend Base;

      using Real = typename TapeTypes::Real;
      using Gradient = typename TapeTypes::Gradient;
      using IndexManager = typename TapeTypes::IndexManager;
      using Identifier = typename TapeTypes::Identifier;

      static_assert(IndexManager::IsLinear, "This class requires an index manager with a linear scheme.");

      JacobianLinearTape() : Base() {}

    protected:

      CODI_INLINE void pushStmtData(Identifier const& index, Config::ArgumentSize const& numberOfArguments) {
        CODI_UNUSED(index);

        this->statementVector.pushData(numberOfArguments);
      }

      CODI_INLINE static void internalEvaluateRevere(
          /* data from call */
          Gradient* adjointVector,
          /* data from jacobian vector */
          size_t& curJacobianPos, size_t const& endJacobianPos, Real const* const rhsJacobians, Identifier const* const rhsIdentifiers ,
          /* data from statement vector */
          size_t& curStmtPos, size_t const& endStmtPos, Config::ArgumentSize const* const numberOfJacobians,
          /* data from index handler */
          size_t const& startAdjointPos, size_t const& endAdjointPos) {

        CODI_UNUSED(endJacobianPos, endStmtPos);

        size_t curAdjointPos = startAdjointPos;

        while(curAdjointPos > endAdjointPos) {

          curStmtPos -= 1;
          Config::ArgumentSize const argsSize = numberOfJacobians[curStmtPos];

          Gradient const lhsAdjoint = adjointVector[curAdjointPos]; // Adjoint positions are shifted since we do not use the zero index

          if(Config::StatementInputTag != argsSize) {
            // No input value, perform regular statement evaluation
            adjointVector[curAdjointPos] = Gradient();

            Base::incrementAdjoints(adjointVector, lhsAdjoint, argsSize, curJacobianPos, rhsJacobians, rhsIdentifiers);
          }

          curAdjointPos -= 1;
        }
      }
  };
}

