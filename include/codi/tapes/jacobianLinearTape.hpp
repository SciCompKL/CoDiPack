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
#include "indices/linearIndexManager.hpp"
#include "interfaces/reverseTapeInterface.hpp"
#include "jacobianBaseTape.hpp"


/** \copydoc codi::Namespace */
namespace codi {


  template<typename _TapeTypes>
  struct JacobianLinearTape : public JacobianBaseTape<_TapeTypes, JacobianLinearTape<_TapeTypes>> {
    public:

      using TapeTypes = DECLARE_DEFAULT(_TapeTypes, TEMPLATE(JacobianTapeTypes<double, double, IndexManagerInterface<int>, ChunkVector>));

      using Base = JacobianBaseTape<_TapeTypes, JacobianLinearTape>;
      friend Base;

      using Real = typename TapeTypes::Real;
      using Gradient = typename TapeTypes::Gradient;
      using IndexManager = typename TapeTypes::IndexManager;
      using Identifier = typename TapeTypes::Identifier;
      using Position = typename Base::Position;

      static_assert(IndexManager::IsLinear, "This class requires an index manager with a linear scheme.");

      JacobianLinearTape() : Base() {}

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

      CODI_INLINE void pushStmtData(Identifier const& index, Config::ArgumentSize const& numberOfArguments) {
        CODI_UNUSED(index);

        this->statementVector.pushData(numberOfArguments);
      }

      template<typename Adjoint>
      CODI_INLINE static void internalEvaluateForward(
          /* data from call */
          Adjoint* adjointVector,
          /* data from jacobian vector */
          size_t& curJacobianPos, size_t const& endJacobianPos, Real const* const rhsJacobians, Identifier const* const rhsIdentifiers ,
          /* data from statement vector */
          size_t& curStmtPos, size_t const& endStmtPos, Config::ArgumentSize const* const numberOfJacobians,
          /* data from index handler */
          size_t const& startAdjointPos, size_t const& endAdjointPos) {

        CODI_UNUSED(endJacobianPos, endStmtPos);

        size_t curAdjointPos = startAdjointPos;

        while(curAdjointPos < endAdjointPos) {

          curAdjointPos += 1;

          Config::ArgumentSize const argsSize = numberOfJacobians[curStmtPos];


          if(Config::StatementInputTag != argsSize) {
            Adjoint lhsAdjoint = Adjoint();

            Base::incrementTangents(adjointVector, lhsAdjoint, argsSize, curJacobianPos, rhsJacobians, rhsIdentifiers);
            adjointVector[curAdjointPos] = lhsAdjoint;
          }


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
          size_t& curStmtPos, size_t const& endStmtPos, Config::ArgumentSize const* const numberOfJacobians,
          /* data from index handler */
          size_t const& startAdjointPos, size_t const& endAdjointPos) {

        CODI_UNUSED(endJacobianPos, endStmtPos);

        size_t curAdjointPos = startAdjointPos;

        while(curAdjointPos > endAdjointPos) {

          curStmtPos -= 1;
          Config::ArgumentSize const argsSize = numberOfJacobians[curStmtPos];

          Adjoint const lhsAdjoint = adjointVector[curAdjointPos]; // Adjoint positions are shifted since we do not use the zero index

          if(Config::StatementInputTag != argsSize) {
            // No input value, perform regular statement evaluation
            adjointVector[curAdjointPos] = Adjoint();

            Base::incrementAdjoints(adjointVector, lhsAdjoint, argsSize, curJacobianPos, rhsJacobians, rhsIdentifiers);
          }

          curAdjointPos -= 1;
        }
      }
  };
}

