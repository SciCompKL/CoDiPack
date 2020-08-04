#pragma once

#include <algorithm>
#include <type_traits>

#include "../aux/macros.hpp"
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

  template<typename _ImplTapeTypes>
  struct JacobianReuseTape : public JacobianBaseTape<_ImplTapeTypes, JacobianReuseTape<_ImplTapeTypes>> {
    public:

      using ImplTapeTypes = CODI_DECLARE_DEFAULT(_ImplTapeTypes, CODI_TEMPLATE(JacobianTapeTypes<double, double, IndexManagerInterface<int>, DefaultChunkedData>));

      using Base = JacobianBaseTape<ImplTapeTypes, JacobianReuseTape>;
      friend Base;

      using Real = typename ImplTapeTypes::Real;
      using Gradient = typename ImplTapeTypes::Gradient;
      using IndexManager = typename ImplTapeTypes::IndexManager;
      using Identifier = typename ImplTapeTypes::Identifier;
      using Position = typename Base::Position;
      using StatementData = typename ImplTapeTypes::StatementData;

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

        using StmtPosition = typename StatementData::Position;
        StmtPosition startStmt = this->externalFunctionData.template extractPosition<StmtPosition>(start);
        StmtPosition endStmt = this->externalFunctionData.template extractPosition<StmtPosition>(end);

        this->statementData.forEachReverse(startStmt, endStmt, clearFunc);
      }

    protected:

      CODI_INLINE void pushStmtData(Identifier const& index, Config::ArgumentSize const& numberOfArguments) {
        this->statementData.pushData(index, numberOfArguments);
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
          /* data from jacobianData */
          size_t& curJacobianPos, size_t const& endJacobianPos, Real const* const rhsJacobians, Identifier const* const rhsIdentifiers ,
          /* data from statementData */
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

