/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2024 Chair for Scientific Computing (SciComp), University of Kaiserslautern-Landau
 * Homepage: http://scicomp.rptu.de
 * Contact:  Prof. Nicolas R. Gauger (codi@scicomp.uni-kl.de)
 *
 * Lead developers: Max Sagebaum, Johannes Blühdorn (SciComp, University of Kaiserslautern-Landau)
 *
 * This file is part of CoDiPack (http://scicomp.rptu.de/software/codi).
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
 *  - SciComp, University of Kaiserslautern-Landau:
 *    - Max Sagebaum
 *    - Johannes Blühdorn
 *    - Former members:
 *      - Tim Albring
 */
#pragma once

#include <algorithm>
#include <type_traits>

#include "../config.h"
#include "../expressions/lhsExpressionInterface.hpp"
#include "../expressions/logic/compileTimeTraversalLogic.hpp"
#include "../expressions/logic/traversalLogic.hpp"
#include "../misc/macros.hpp"
#include "../traits/expressionTraits.hpp"
#include "data/chunk.hpp"
#include "indices/indexManagerInterface.hpp"
#include "interfaces/editingTapeInterface.hpp"
#include "interfaces/reverseTapeInterface.hpp"
#include "jacobianBaseTape.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Final implementation for a Jacobian tape with a reuse index management.
   *
   * This class implements the interface methods from the JacobianBaseTape.
   *
   * @tparam T_TapeTypes  JacobianTapeTypes definition.
   */
  template<typename T_TapeTypes>
  struct JacobianReuseTape
      : public JacobianBaseTape<T_TapeTypes, JacobianReuseTape<T_TapeTypes>>,
        public EditingTapeInterface<typename JacobianBaseTape<T_TapeTypes, JacobianReuseTape<T_TapeTypes>>::Position> {
    public:

      using TapeTypes = CODI_DD(T_TapeTypes,
                                CODI_T(JacobianTapeTypes<double, double, IndexManagerInterface<int>,
                                                         DefaultChunkedData>));  ///< See JacobianReuseTape.

      using Base = JacobianBaseTape<T_TapeTypes, JacobianReuseTape>;  ///< Base class abbreviation.
      friend Base;  ///< Allow the base class to call protected and private methods.

      using Real = typename TapeTypes::Real;                    ///< See TapeTypesInterface.
      using Gradient = typename TapeTypes::Gradient;            ///< See TapeTypesInterface.
      using IndexManager = typename TapeTypes::IndexManager;    ///< See TapeTypesInterface.
      using Identifier = typename TapeTypes::Identifier;        ///< See TapeTypesInterface.
      using Position = typename Base::Position;                 ///< See TapeTypesInterface.
      using StatementData = typename TapeTypes::StatementData;  ///< See JacobianTapeTypes.

      CODI_STATIC_ASSERT(!IndexManager::IsLinear, "This class requires an index manager with a reuse scheme.");

      /// Constructor
      JacobianReuseTape() : Base() {}

      /*******************************************************************************/
      /// @name Missing functions from FullTapeInterface
      /// @{

      using Base::clearAdjoints;

      /// \copydoc codi::PositionalEvaluationTapeInterface::clearAdjoints
      void clearAdjoints(Position const& start, Position const& end,
                         AdjointsManagement adjointsManagement = AdjointsManagement::Automatic) {
        if (AdjointsManagement::Automatic == adjointsManagement) {
          this->adjoints.beginUse();
        }

        Identifier adjointsSize = (Identifier)this->adjoints.size();

        auto clearFunc = [&adjointsSize, this](Identifier* index, Config::ArgumentSize* stmtSize) {
          CODI_UNUSED(stmtSize);

          if (*index < adjointsSize) {
            this->adjoints[*index] = Gradient();
          }
        };

        using StmtPosition = typename StatementData::Position;
        StmtPosition startStmt = this->llfByteData.template extractPosition<StmtPosition>(start);
        StmtPosition endStmt = this->llfByteData.template extractPosition<StmtPosition>(end);

        this->statementData.forEachReverse(startStmt, endStmt, clearFunc);

        if (AdjointsManagement::Automatic == adjointsManagement) {
          this->adjoints.endUse();
        }
      }

      /// @}

    protected:

      /*******************************************************************************/
      /// @name Functions from JacobianBaseTape
      /// @{

      /// \copydoc codi::JacobianBaseTape::pushStmtData <br><br>
      /// Both arguments are pushed to the tape.
      CODI_INLINE void pushStmtData(Identifier const& index, Config::ArgumentSize const& numberOfArguments) {
        this->statementData.pushData(index, numberOfArguments);
      }

      /// \copydoc codi::JacobianBaseTape::internalEvaluateForward_EvalStatements
      template<typename Adjoint>
      CODI_INLINE static void internalEvaluateForward_EvalStatements(
          /* data from call */
          JacobianReuseTape& tape, Adjoint* adjointVector,
          /* data from low level function byte data vector */
          size_t& curLLFByteDataPos, size_t const& endLLFByteDataPos, char* dataPtr,
          /* data from low level function info data vector */
          size_t& curLLFInfoDataPos, size_t const& endLLFInfoDataPos, Config::LowLevelFunctionToken* const tokenPtr,
          Config::LowLevelFunctionDataSize* const dataSizePtr,
          /* data from jacobian vector */
          size_t& curJacobianPos, size_t const& endJacobianPos, Real const* const rhsJacobians,
          Identifier const* const rhsIdentifiers,
          /* data from statement vector */
          size_t& curStmtPos, size_t const& endStmtPos, Identifier const* const lhsIdentifiers,
          Config::ArgumentSize const* const numberOfJacobians) {
        CODI_UNUSED(endJacobianPos, endLLFByteDataPos, endLLFInfoDataPos);

        typename Base::template VectorAccess<Adjoint> vectorAccess(adjointVector);

        while (curStmtPos < endStmtPos) CODI_Likely {
          Config::ArgumentSize const argsSize = numberOfJacobians[curStmtPos];

          if (Config::StatementLowLevelFunctionTag == argsSize) CODI_Unlikely {
            Base::template callLowLevelFunction<LowLevelFunctionEntryCallKind::Forward>(
                tape, true, curLLFByteDataPos, dataPtr, curLLFInfoDataPos, tokenPtr, dataSizePtr, &vectorAccess);
          } else CODI_Likely {
            Adjoint lhsAdjoint = Adjoint();
            Base::incrementTangents(adjointVector, lhsAdjoint, argsSize, curJacobianPos, rhsJacobians, rhsIdentifiers);

            adjointVector[lhsIdentifiers[curStmtPos]] = lhsAdjoint;

            EventSystem<JacobianReuseTape>::notifyStatementEvaluateListeners(
                tape, lhsIdentifiers[curStmtPos], GradientTraits::dim<Adjoint>(),
                GradientTraits::toArray(lhsAdjoint).data());
          }

          curStmtPos += 1;
        }
      }

      /// \copydoc codi::JacobianBaseTape::internalEvaluateReverse_EvalStatements
      template<typename Adjoint>
      CODI_INLINE static void internalEvaluateReverse_EvalStatements(
          /* data from call */
          JacobianReuseTape& tape, Adjoint* adjointVector,
          /* data from low level function byte data vector */
          size_t& curLLFByteDataPos, size_t const& endLLFByteDataPos, char* dataPtr,
          /* data from low level function info data vector */
          size_t& curLLFInfoDataPos, size_t const& endLLFInfoDataPos, Config::LowLevelFunctionToken* const tokenPtr,
          Config::LowLevelFunctionDataSize* const dataSizePtr,
          /* data from jacobianData */
          size_t& curJacobianPos, size_t const& endJacobianPos, Real const* const rhsJacobians,
          Identifier const* const rhsIdentifiers,
          /* data from statementData */
          size_t& curStmtPos, size_t const& endStmtPos, Identifier const* const lhsIdentifiers,
          Config::ArgumentSize const* const numberOfJacobians) {
        CODI_UNUSED(endJacobianPos, endLLFByteDataPos, endLLFInfoDataPos);

        typename Base::template VectorAccess<Adjoint> vectorAccess(adjointVector);

        while (curStmtPos > endStmtPos) CODI_Likely {
          curStmtPos -= 1;

          Config::ArgumentSize const argsSize = numberOfJacobians[curStmtPos];

          if (Config::StatementLowLevelFunctionTag == argsSize) CODI_Unlikely {
            Base::template callLowLevelFunction<LowLevelFunctionEntryCallKind::Reverse>(
                tape, false, curLLFByteDataPos, dataPtr, curLLFInfoDataPos, tokenPtr, dataSizePtr, &vectorAccess);
          } else CODI_Likely {
            Adjoint const lhsAdjoint = adjointVector[lhsIdentifiers[curStmtPos]];

            EventSystem<JacobianReuseTape>::notifyStatementEvaluateListeners(
                tape, lhsIdentifiers[curStmtPos], GradientTraits::dim<Adjoint>(),
                GradientTraits::toArray(lhsAdjoint).data());

            adjointVector[lhsIdentifiers[curStmtPos]] = Adjoint();
            Base::incrementAdjoints(adjointVector, lhsAdjoint, argsSize, curJacobianPos, rhsJacobians, rhsIdentifiers);
          }
        }
      }

      /// @}

    public:

      /*******************************************************************************/
      /// @name Functions from EditingTapeInterface
      /// @{

      /// \copydoc codi::EditingTapeInterface::erase(T_Position const& start, T_Position const& end) <br>
      /// Implementation: Instantiates a temporary tape. If called often, this can become a bottleneck. The variant of
      /// erase that takes a reference to a helper tape should be used in this case.
      CODI_INLINE void erase(Position const& start, Position const& end) {
        JacobianReuseTape emptyTape;
        erase(start, end, emptyTape);
      }

      // clang-format off
      /// \copydoc codi::EditingTapeInterface::erase(T_Position const& start, T_Position const& end, EditingTapeInterface& emptyTape)
      // clang-format on
      CODI_INLINE void erase(Position const& start, Position const& end, JacobianReuseTape& emptyTape) {
        // Store the tail after the part to be erased in the helper tape.
        emptyTape.append(*this, end, this->getPosition());
        // Reset the tape to before the erased part and re-append the tail. This accounts for external function position
        // correction.
        this->resetTo(start);
        this->append(emptyTape, emptyTape.getZeroPosition(), emptyTape.getPosition());
        emptyTape.reset();
      }

      /// \copydoc codi::EditingTapeInterface::append
      CODI_INLINE void append(JacobianReuseTape& srcTape, Position const& start, Position const& end) {
        srcTape.llfByteData.evaluateForward(start, end, JacobianReuseTape::internalAppend, this);
      }

      /// @}

    private:

      static CODI_INLINE void internalAppend(
          /* data from call */
          JacobianReuseTape* dstTape,
          /* data from low level function byte data vector */
          size_t& curLLFByteDataPos, size_t const& endLLFByteDataPos, char* dataPtr,
          /* data from low level function info data vector */
          size_t& curLLFInfoDataPos, size_t const& endLLFInfoDataPos, Config::LowLevelFunctionToken* const tokenPtr,
          Config::LowLevelFunctionDataSize* const dataSizePtr,
          /* data from jacobianData */
          size_t& curJacobianPos, size_t const& endJacobianPos, Real const* const rhsJacobians,
          Identifier const* const rhsIdentifiers,
          /* data from statementData */
          size_t& curStmtPos, size_t const& endStmtPos, Identifier const* const lhsIdentifiers,
          Config::ArgumentSize const* const numberOfJacobians) {
        CODI_UNUSED(endLLFByteDataPos, endLLFInfoDataPos, endJacobianPos);

        while (curStmtPos < endStmtPos) {
          Config::ArgumentSize const argsSize = numberOfJacobians[curStmtPos];
          if (Config::StatementLowLevelFunctionTag == argsSize) CODI_Unlikely {
            Config::LowLevelFunctionToken token = tokenPtr[curLLFInfoDataPos];
            size_t dataSize = dataSizePtr[curLLFInfoDataPos];

            // Create the store on the new tape.
            ByteDataView dstDataStore = {};
            dstTape->pushLowLevelFunction(token, dataSize, dstDataStore);

            // Copy the data.
            dstDataStore.write(&dataPtr[curLLFByteDataPos], dataSize);

            curLLFInfoDataPos += 1;
            curLLFByteDataPos += dataSize;
          } else CODI_Likely {
            // Manual statement push.
            dstTape->statementData.reserveItems(1);
            dstTape->jacobianData.reserveItems(numberOfJacobians[curStmtPos]);

            dstTape->pushStmtData(lhsIdentifiers[curStmtPos], numberOfJacobians[curStmtPos]);
            size_t curJacobianEnd = curJacobianPos + numberOfJacobians[curStmtPos];

            while (curJacobianPos < curJacobianEnd) {
              dstTape->jacobianData.pushData(rhsJacobians[curJacobianPos], rhsIdentifiers[curJacobianPos]);
              ++curJacobianPos;
            }
          }

          ++curStmtPos;
        }
      }
  };
}
