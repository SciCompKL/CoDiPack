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
#include "indices/linearIndexManager.hpp"
#include "interfaces/reverseTapeInterface.hpp"
#include "jacobianBaseTape.hpp"


/** \copydoc codi::Namespace */
namespace codi {


  /**
   * @brief Final implementation for a Jacobian tape with a linear index management.
   *
   * This class implements the interface methods from the JacobianBaseTape.
   *
   * @tparam _TapeTypes  JacobianTapeTypes definition.
   */
  template<typename _TapeTypes>
  struct JacobianLinearTape : public JacobianBaseTape<_TapeTypes, JacobianLinearTape<_TapeTypes>> {
    public:

      using TapeTypes = CODI_DD(_TapeTypes, CODI_T(JacobianTapeTypes<double, double,
                        IndexManagerInterface<int>, DefaultChunkedData>)); ///< See JacobianLinearTape

      using Base = JacobianBaseTape<TapeTypes, JacobianLinearTape>; ///< Base class abbreviation
      friend Base; ///< Allow the base class to call protected and private methods.

      using Real = typename TapeTypes::Real;                    ///< See TapeTypesInterface.
      using Gradient = typename TapeTypes::Gradient;            ///< See TapeTypesInterface.
      using IndexManager = typename TapeTypes::IndexManager;    ///< See TapeTypesInterface.
      using Identifier = typename TapeTypes::Identifier;        ///< See TapeTypesInterface.
      using Position = typename Base::Position;                 ///< See TapeTypesInterface.

      static_assert(IndexManager::IsLinear, "This class requires an index manager with a linear scheme.");

      /// Constructor
      JacobianLinearTape() : Base() {}

      using Base::clearAdjoints;

      /// \copydoc codi::PositionalEvaluationTapeInterface::clearAdjoints
      void clearAdjoints(Position const& start, Position const& end) {

        using IndexPosition = typename IndexManager::Position;
        IndexPosition startIndex = this->externalFunctionData.template extractPosition<IndexPosition>(start);
        IndexPosition endIndex = this->externalFunctionData.template extractPosition<IndexPosition>(end);

        startIndex = std::min(startIndex, (IndexPosition)this->adjoints.size() - 1);
        endIndex = std::min(endIndex, (IndexPosition)this->adjoints.size() - 1);

        for(IndexPosition curPos = endIndex + 1; curPos <= startIndex; curPos += 1) {
          this->adjoints[curPos] = Gradient();
        }
      }

    protected:

      /// Only number of arguments is required for linear index managers
      CODI_INLINE void pushStmtData(Identifier const& index, Config::ArgumentSize const& numberOfArguments) {
        CODI_UNUSED(index);

        this->statementData.pushData(numberOfArguments);
      }

      /// \copydoc codi::JacobianBaseTape::internalEvaluateForwardStack
      template<typename Adjoint>
      CODI_INLINE static void internalEvaluateForwardStack(
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

      /// \copydoc codi::JacobianBaseTape::internalEvaluateReverseStack
      template<typename Adjoint>
      CODI_INLINE static void internalEvaluateReverseStack(
          /* data from call */
          Adjoint* adjointVector,
          /* data from jacobianData */
          size_t& curJacobianPos, size_t const& endJacobianPos, Real const* const rhsJacobians, Identifier const* const rhsIdentifiers ,
          /* data from statementData */
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

