/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2021 Chair for Scientific Computing (SciComp), TU Kaiserslautern
 * Homepage: http://www.scicomp.uni-kl.de
 * Contact:  Prof. Nicolas R. Gauger (codi@scicomp.uni-kl.de)
 *
 * Lead developers: Max Sagebaum, Johannes Blühdorn (SciComp, TU Kaiserslautern)
 *
 * This file is part of CoDiPack (http://www.scicomp.uni-kl.de/software/codi).
 *
 * CoDiPack is free software: you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * CoDiPack is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty
 * of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 * See the GNU General Public License for more details.
 * You should have received a copy of the GNU
 * General Public License along with CoDiPack.
 * If not, see <http://www.gnu.org/licenses/>.
 *
 * For other licensing options please contact us.
 *
 * Authors:
 *  - SciComp, TU Kaiserslautern:
 *    - Max Sagebaum
 *    - Johannes Blühdorn
 *    - Former members:
 *      - Tim Albring
 */
#pragma once

#include <algorithm>
#include <type_traits>

#include "../misc/macros.hpp"
#include "../config.h"
#include "../expressions/lhsExpressionInterface.hpp"
#include "../expressions/logic/compileTimeTraversalLogic.hpp"
#include "../expressions/logic/traversalLogic.hpp"
#include "../traits/expressionTraits.hpp"
#include "data/chunk.hpp"
#include "indices/indexManagerInterface.hpp"
#include "interfaces/reverseTapeInterface.hpp"
#include "parallelJacobianBaseTape.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Final implementation for a parallel Jacobian tape with a reuse index management.
   *
   * This class implements the interface methods from the ParallelJacobianBaseTape.
   *
   * @tparam T_TapeTypes  JacobianTapeTypes definition.
   */
  template<typename T_TapeTypes>
  struct ParallelJacobianReuseTape : public ParallelJacobianBaseTape<T_TapeTypes,
                                                                     ParallelJacobianReuseTape<T_TapeTypes>> {
    public:

      using TapeTypes = CODI_DD(T_TapeTypes,
                                CODI_T(JacobianTapeTypes<double, double, IndexManagerInterface<int>,
                                                         DefaultChunkedData>));  ///< See JacobianReuseTape.

      using Base = ParallelJacobianBaseTape<TapeTypes, ParallelJacobianReuseTape>;  ///< Base class abbreviation.
      friend Base;  ///< Allow the base class to call protected and private methods.
      friend typename Base::Base;  ///< Allow the base class to call protected and private methods.

      using Real = typename TapeTypes::Real;                    ///< See TapeTypesInterface.
      using Gradient = typename TapeTypes::Gradient;            ///< See TapeTypesInterface.
      using IndexManager = typename TapeTypes::IndexManager;    ///< See TapeTypesInterface.
      using Identifier = typename TapeTypes::Identifier;        ///< See TapeTypesInterface.
      using Position = typename Base::Position;                 ///< See TapeTypesInterface.
      using StatementData = typename TapeTypes::StatementData;  ///< See JacobianTapeTypes.

      static_assert(!IndexManager::IsLinear, "This class requires an index manager with a reuse scheme.");

      /// Constructor
      ParallelJacobianReuseTape() : Base() {}

      /// \copydoc codi::ReverseTapeInterface::getTapeValues()
      TapeValues internalGetTapeValues() {
        TapeValues values = TapeValues("CoDi Tape Statistics ( ParallelJacobianReuseTape )");
        Base::internalAddTapeValues(values);
        return values;
      }

      using Base::clearAdjoints;

      /// \copydoc codi::PositionalEvaluationTapeInterface::clearAdjoints
      void clearAdjoints(Position const& start, Position const& end) {
        auto clearFunc = [this](Identifier* index, Config::ArgumentSize* stmtSize) {
          CODI_UNUSED(stmtSize);

          if (*index < (Identifier)this->adjoints.size()) {
            this->adjoints[*index] = Gradient();
          }
        };

        Base::LockForUse lock(Base::adjointsMutex);

        using StmtPosition = typename StatementData::Position;
        StmtPosition startStmt = this->externalFunctionData.template extractPosition<StmtPosition>(start);
        StmtPosition endStmt = this->externalFunctionData.template extractPosition<StmtPosition>(end);

        this->statementData.forEachReverse(startStmt, endStmt, clearFunc);
      }

    protected:

      /// \copydoc codi::JacobianBaseTape::pushStmtData
      /// Both arguments are pushed to the tape.
      CODI_INLINE void pushStmtData(Identifier const& index, Config::ArgumentSize const& numberOfArguments) {
        this->statementData.pushData(index, numberOfArguments);
      }

      /// \copydoc codi::JacobianBaseTape::internalEvaluateForward_Step3_EvalStatements
      template<typename Adjoint>
      CODI_INLINE static void internalEvaluateForward_Step3_EvalStatements(
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

      /// \copydoc codi::JacobianBaseTape::internalEvaluateReverse_Step3_EvalStatements
      template<typename Adjoint>
      CODI_INLINE static void internalEvaluateReverse_Step3_EvalStatements(
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
