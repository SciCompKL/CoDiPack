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
#include "indices/indexManagerInterface.hpp"
#include "interfaces/reverseTapeInterface.hpp"
#include "jacobianBaseTape.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  template<typename _TapeTypes>
  struct JacobianReuseTape : public JacobianBaseTape<_TapeTypes, JacobianReuseTape<_TapeTypes>> {
    public:

      using TapeTypes = DECLARE_DEFAULT(_TapeTypes, TEMPLATE(JacobianTapeTypes<double, double, IndexManagerInterface<int>));

      using Base = JacobianBaseTape<_TapeTypes, JacobianReuseTape>;
      friend Base;

      using Real = typename TapeTypes::Real;
      using Gradient = typename TapeTypes::Gradient;
      using IndexManager = typename TapeTypes::IndexManager;
      using Identifier = typename TapeTypes::Identifier;

      static_assert(!IndexManager::IsLinear, "This class requires an index manager with a reuse scheme.");

      JacobianReuseTape() : Base() {}

      template<typename Lhs>
      CODI_INLINE void registerInput(LhsExpressionInterface<Real, Gradient, JacobianReuseTape, Lhs>& value) {
        Base::indexManager.get().assignUnusedIndex(value.cast().getIdentifier());
      }

    protected:

      CODI_INLINE void pushStmtData(Identifier const& index, Config::ArgumentSize const& numberOfArguments) {
        this->statementVector.pushData(index, numberOfArguments);
      }

      CODI_INLINE static void internalEvaluateRevere(
          /* data from call */
          Gradient* adjointVector,
          /* data from jacobian vector */
          size_t& curJacobianPos, size_t const& endJacobianPos, Real const* const rhsJacobians, Identifier const* const rhsIdentifiers ,
          /* data from statement vector */
          size_t& curStmtPos, size_t const& endStmtPos, Identifier const* const lhsIdentifiers, Config::ArgumentSize const* const numberOfJacobians) {

        CODI_UNUSED(endJacobianPos);

        while(curStmtPos > endStmtPos) {
          curStmtPos -= 1;

          Gradient const lhsAdjoint = adjointVector[lhsIdentifiers[curStmtPos]];
          adjointVector[lhsIdentifiers[curStmtPos]] = Gradient();

          Base::incrementAdjoints(adjointVector, lhsAdjoint, numberOfJacobians[curStmtPos], curJacobianPos, rhsJacobians, rhsIdentifiers);
        }
      }
  };
}

