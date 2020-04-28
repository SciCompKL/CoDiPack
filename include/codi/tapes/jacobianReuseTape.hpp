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
#include "indices/indexManagerInterface.hpp"
#include "interfaces/reverseTapeInterface.hpp"
#include "jacobianBaseTape.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  template<typename _TapeTypes>
  struct JacobianReuseTape : public JacobianBaseTape<_TapeTypes, JacobianReuseTape<_TapeTypes>> {
    public:

      using TapeTypes = DECLARE_DEFAULT(_TapeTypes, TEMPLATE(JacobianTapeTypes<double, double, IndexManagerInterface<int>, ChunkVector>));

      using Base = JacobianBaseTape<_TapeTypes, JacobianReuseTape>;
      friend Base;

      using Real = typename TapeTypes::Real;
      using Gradient = typename TapeTypes::Gradient;
      using IndexManager = typename TapeTypes::IndexManager;
      using Identifier = typename TapeTypes::Identifier;
      using Position = typename Base::Position;
      using StatementVector = typename TapeTypes::StatementVector;

      static_assert(!IndexManager::IsLinear, "This class requires an index manager with a reuse scheme.");

      JacobianReuseTape() : Base() {}

      using Base::clearAdjoints;
      void clearAdjoints(Position const& start, Position const& end) {

        // clear adjoints
        auto clearFunc = [this] (Config::ArgumentSize* stmtSize, Identifier* index) {
          CODI_UNUSED(stmtSize);

          if(*index < this->adjoints.size()) {
            this->adjoints[*index] = Gradient();
          }
        };

        using StmtPosition = typename StatementVector::Position;
        StmtPosition startStmt = this->externalFunctionVector.template extractPosition<StmtPosition>(start);
        StmtPosition endStmt = this->externalFunctionVector.template extractPosition<StmtPosition>(end);

        this->statementVector.forEachReverse(startStmt, endStmt, clearFunc);
      }

    protected:

      CODI_INLINE void pushStmtData(Identifier const& index, Config::ArgumentSize const& numberOfArguments) {
        this->statementVector.pushData(index, numberOfArguments);
      }

      template<typename Adjoint>
      CODI_INLINE static void internalEvaluateForward(
          /* data from call */
          Adjoint* adjointVector,
          /* data from jacobian vector */
          size_t& curJacobianPos, size_t const& endJacobianPos, Real const* const rhsJacobians, Identifier const* const rhsIdentifiers ,
          /* data from statement vector */
          size_t& curStmtPos, size_t const& endStmtPos, Identifier const* const lhsIdentifiers, Config::ArgumentSize const* const numberOfJacobians) {

        CODI_UNUSED(endJacobianPos);

        while(curStmtPos < endStmtPos) {

          Adjoint lhsAdjoint = Adjoint();
          Base::incrementTangents(adjointVector, lhsAdjoint, numberOfJacobians[curStmtPos], curJacobianPos, rhsJacobians, rhsIdentifiers);

          adjointVector[lhsIdentifiers[curStmtPos]] = lhsAdjoint;

          curStmtPos += 1;
        }
      }

      template<typename Adjoint>
      CODI_INLINE static void internalEvaluateReverse(
          /* data from call */
          Adjoint* adjointVector,
          /* data from jacobian vector */
          size_t& curJacobianPos, size_t const& endJacobianPos, Real const* const rhsJacobians, Identifier const* const rhsIdentifiers ,
          /* data from statement vector */
          size_t& curStmtPos, size_t const& endStmtPos, Identifier const* const lhsIdentifiers, Config::ArgumentSize const* const numberOfJacobians) {

        CODI_UNUSED(endJacobianPos);

        while(curStmtPos > endStmtPos) {
          curStmtPos -= 1;

          Adjoint const lhsAdjoint = adjointVector[lhsIdentifiers[curStmtPos]];
          adjointVector[lhsIdentifiers[curStmtPos]] = Adjoint();

          Base::incrementAdjoints(adjointVector, lhsAdjoint, numberOfJacobians[curStmtPos], curJacobianPos, rhsJacobians, rhsIdentifiers);
        }
      }
  };
}

