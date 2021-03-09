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

  /**
   * @brief Final implementation for a Jacobian tape with a reuse index management.
   *
   * This class implements the interface methods from the JacobianBaseTape.
   *
   * @tparam _TapeTypes  JacobianTapeTypes definition.
   */
  template<typename _TapeTypes>
  struct JacobianReuseTape : public JacobianBaseTape<_TapeTypes, JacobianReuseTape<_TapeTypes>> {
    public:

      using TapeTypes = CODI_DD(_TapeTypes, CODI_T(JacobianTapeTypes<double, double, IndexManagerInterface<int>,
                                                                     DefaultChunkedData>));  ///< See JacobianReuseTape

      using Base = JacobianBaseTape<TapeTypes, JacobianReuseTape>;  ///< Base class abbreviation
      friend Base;  ///< Allow the base class to call protected and private methods.

      using Real = typename TapeTypes::Real;                    ///< See TapeTypesInterface.
      using Gradient = typename TapeTypes::Gradient;            ///< See TapeTypesInterface.
      using IndexManager = typename TapeTypes::IndexManager;    ///< See TapeTypesInterface.
      using Identifier = typename TapeTypes::Identifier;        ///< See TapeTypesInterface.
      using Position = typename Base::Position;                 ///< See TapeTypesInterface.
      using StatementData = typename TapeTypes::StatementData;  ///< See JacobianTapeTypes

      static_assert(!IndexManager::IsLinear, "This class requires an index manager with a reuse scheme.");

      /// Constructor
      JacobianReuseTape() : Base() {}

      using Base::clearAdjoints;

      /// \copydoc codi::PositionalEvaluationTapeInterface::clearAdjoints
      void clearAdjoints(Position const& start, Position const& end) {
        // clear adjoints
        auto clearFunc = [this](Identifier* index, Config::ArgumentSize* stmtSize) {
          CODI_UNUSED(stmtSize);

          if (*index < (Identifier)this->adjoints.size()) {
            this->adjoints[*index] = Gradient();
          }
        };

        using StmtPosition = typename StatementData::Position;
        StmtPosition startStmt = this->externalFunctionData.template extractPosition<StmtPosition>(start);
        StmtPosition endStmt = this->externalFunctionData.template extractPosition<StmtPosition>(end);

        this->statementData.forEachReverse(startStmt, endStmt, clearFunc);
      }

    protected:

      /// Both arguments are pushed to the tape.
      CODI_INLINE void pushStmtData(Identifier const& index, Config::ArgumentSize const& numberOfArguments) {
        this->statementData.pushData(index, numberOfArguments);
      }

      /// \copydoc codi::JacobianBaseTape::internalEvaluateForwardStack
      template<typename Adjoint>
      CODI_INLINE static void internalEvaluateForwardStack(
          /* data from call */
          Adjoint* adjointVector,
          /* data from jacobian vector */
          size_t& curJacobianPos, size_t const& endJacobianPos, Real const* const rhsJacobians,
          Identifier const* const rhsIdentifiers,
          /* data from statement vector */
          size_t& curStmtPos, size_t const& endStmtPos, Identifier const* const lhsIdentifiers,
          Config::ArgumentSize const* const numberOfJacobians) {
        CODI_UNUSED(endJacobianPos);

        while (curStmtPos < endStmtPos) {
          Adjoint lhsAdjoint = Adjoint();
          Base::incrementTangents(adjointVector, lhsAdjoint, numberOfJacobians[curStmtPos], curJacobianPos,
                                  rhsJacobians, rhsIdentifiers);

          adjointVector[lhsIdentifiers[curStmtPos]] = lhsAdjoint;

          curStmtPos += 1;
        }
      }

      /// \copydoc codi::JacobianBaseTape::internalEvaluateReverseStack
      template<typename Adjoint>
      CODI_INLINE static void internalEvaluateReverseStack(
          /* data from call */
          Adjoint* adjointVector,
          /* data from jacobianData */
          size_t& curJacobianPos, size_t const& endJacobianPos, Real const* const rhsJacobians,
          Identifier const* const rhsIdentifiers,
          /* data from statementData */
          size_t& curStmtPos, size_t const& endStmtPos, Identifier const* const lhsIdentifiers,
          Config::ArgumentSize const* const numberOfJacobians) {
        CODI_UNUSED(endJacobianPos);

        while (curStmtPos > endStmtPos) {
          curStmtPos -= 1;

          Adjoint const lhsAdjoint = adjointVector[lhsIdentifiers[curStmtPos]];
          adjointVector[lhsIdentifiers[curStmtPos]] = Adjoint();

          Base::incrementAdjoints(adjointVector, lhsAdjoint, numberOfJacobians[curStmtPos], curJacobianPos,
                                  rhsJacobians, rhsIdentifiers);
        }
      }
  };
}
