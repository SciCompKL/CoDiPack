/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2025 Chair for Scientific Computing (SciComp), University of Kaiserslautern-Landau
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
#include "../traits/adjointVectorTraits.hpp"
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
   * @tparam T_TapeTypes  JacobianTapeTypes definition.
   */
  template<typename T_TapeTypes>
  struct JacobianLinearTape : public JacobianBaseTape<T_TapeTypes, JacobianLinearTape<T_TapeTypes>> {
    public:

      using TapeTypes = CODI_DD(T_TapeTypes,
                                CODI_T(JacobianTapeTypes<double, double, IndexManagerInterface<int>,
                                                         DefaultChunkedData>));  ///< See JacobianLinearTape.

      using Base = JacobianBaseTape<TapeTypes, JacobianLinearTape>;  ///< Base class abbreviation.
      friend Base;  ///< Allow the base class to call protected and private methods.

      using Real = typename TapeTypes::Real;                              ///< See TapeTypesInterface.
      using Gradient = typename TapeTypes::Gradient;                      ///< See TapeTypesInterface.
      using IndexManager = typename TapeTypes::IndexManager;              ///< See TapeTypesInterface.
      using Identifier = typename TapeTypes::Identifier;                  ///< See TapeTypesInterface.
      using ActiveTypeTapeData = typename TapeTypes::ActiveTypeTapeData;  ///< See TapeTypesInterface.
      using Position = typename Base::Position;                           ///< See TapeTypesInterface.

      CODI_STATIC_ASSERT(IndexManager::IsLinear, "This class requires an index manager with a linear scheme.");

      /// Constructor
      JacobianLinearTape() : Base() {}

      using Base::clearAdjoints;

      /// \copydoc codi::PositionalEvaluationTapeInterface::clearAdjoints
      void clearAdjoints(Position const& start, Position const& end,
                         AdjointsManagement adjointsManagement = AdjointsManagement::Automatic) {
        if (AdjointsManagement::Automatic == adjointsManagement) {
          this->adjoints.beginUse();
        }

        using IndexPosition = CODI_DD(typename IndexManager::Position, int);
        IndexPosition startIndex = this->llfByteData.template extractPosition<IndexPosition>(start);
        IndexPosition endIndex = this->llfByteData.template extractPosition<IndexPosition>(end);

        startIndex = std::min(startIndex, (IndexPosition)this->adjoints.size() - 1);
        endIndex = std::min(endIndex, (IndexPosition)this->adjoints.size() - 1);

        for (IndexPosition curPos = endIndex + 1; curPos <= startIndex; curPos += 1) {
          this->adjoints[curPos] = Gradient();
        }

        if (AdjointsManagement::Automatic == adjointsManagement) {
          this->adjoints.endUse();
        }
      }

      /// \copydoc codi::CustomAdjointVectorEvaluationTapeInterface::clearCustomAdjoints
      template<typename AdjointVector>
      void clearCustomAdjoints(Position const& start, Position const& end, AdjointVector&& data) {
        using IndexPosition = CODI_DD(typename IndexManager::Position, int);
        IndexPosition startIndex = this->llfByteData.template extractPosition<IndexPosition>(start);
        IndexPosition endIndex = this->llfByteData.template extractPosition<IndexPosition>(end);

        for (IndexPosition curPos = endIndex + 1; curPos <= startIndex; curPos += 1) {
          data[curPos] = AdjointVectorTraits::Gradient<AdjointVector>();
        }
      }

    protected:

      /// \copydoc codi::JacobianBaseTape::pushStmtData <br><br>
      /// Only the number of arguments is required for linear index managers.
      CODI_INLINE void pushStmtData(Identifier const& index, Config::ArgumentSize const& numberOfArguments) {
        CODI_UNUSED(index);

        this->statementData.pushData(numberOfArguments);
      }

      /// \copydoc codi::JacobianBaseTape::internalEvaluateForward_EvalStatements
      template<typename AdjointVector>
      CODI_INLINE static void internalEvaluateForward_EvalStatements(
          /* data from call */
          JacobianLinearTape& tape, AdjointVector&& adjointVector,
          /* data from low level function byte data vector */
          size_t& curLLFByteDataPos, size_t const& endLLFByteDataPos, char* dataPtr,
          /* data from low level function info data vector */
          size_t& curLLFInfoDataPos, size_t const& endLLFInfoDataPos, Config::LowLevelFunctionToken* const tokenPtr,
          Config::LowLevelFunctionDataSize* const dataSizePtr,
          /* data from jacobian vector */
          size_t& curJacobianPos, size_t const& endJacobianPos, Real const* const rhsJacobians,
          Identifier const* const rhsIdentifiers,
          /* data from statement vector */
          size_t& curStmtPos, size_t const& endStmtPos, Config::ArgumentSize const* const numberOfJacobians,
          /* data from index handler */
          size_t const& startAdjointPos, size_t const& endAdjointPos) {
        CODI_UNUSED(endJacobianPos, endStmtPos, endLLFByteDataPos, endLLFInfoDataPos);

        using Adjoint = AdjointVectorTraits::Gradient<AdjointVector>;

        typename Base::template VectorAccess<AdjointVector> vectorAccess(adjointVector);

        size_t curAdjointPos = startAdjointPos;

        while (curAdjointPos < endAdjointPos) CODI_Likely {
          curAdjointPos += 1;

          Config::ArgumentSize const argsSize = numberOfJacobians[curStmtPos];

          if (Config::StatementLowLevelFunctionTag == argsSize) CODI_Unlikely {
            Base::template callLowLevelFunction<LowLevelFunctionEntryCallKind::Forward>(
                tape, true, curLLFByteDataPos, dataPtr, curLLFInfoDataPos, tokenPtr, dataSizePtr, &vectorAccess);
          } else if (Config::StatementInputTag == argsSize) CODI_Unlikely {
            // Do nothing.
          } else CODI_Likely {
            Adjoint lhsAdjoint = Adjoint();
            Base::incrementTangents(adjointVector, lhsAdjoint, argsSize, curJacobianPos, rhsJacobians, rhsIdentifiers);
            adjointVector[curAdjointPos] = lhsAdjoint;

            EventSystem<JacobianLinearTape>::notifyStatementEvaluateListeners(
                tape, curAdjointPos, GradientTraits::dim<Adjoint>(), GradientTraits::toArray(lhsAdjoint).data());
          }

          curStmtPos += 1;
        }
      }

      /// \copydoc codi::JacobianBaseTape::internalEvaluateReverse_EvalStatements
      template<typename AdjointVector>
      CODI_INLINE static void internalEvaluateReverse_EvalStatements(
          /* data from call */
          JacobianLinearTape& tape, AdjointVector&& adjointVector,
          /* data from low level function byte data vector */
          size_t& curLLFByteDataPos, size_t const& endLLFByteDataPos, char* dataPtr,
          /* data from low level function info data vector */
          size_t& curLLFInfoDataPos, size_t const& endLLFInfoDataPos, Config::LowLevelFunctionToken* const tokenPtr,
          Config::LowLevelFunctionDataSize* const dataSizePtr,
          /* data from jacobianData */
          size_t& curJacobianPos, size_t const& endJacobianPos, Real const* const rhsJacobians,
          Identifier const* const rhsIdentifiers,
          /* data from statementData */
          size_t& curStmtPos, size_t const& endStmtPos, Config::ArgumentSize const* const numberOfJacobians,
          /* data from index handler */
          size_t const& startAdjointPos, size_t const& endAdjointPos) {
        CODI_UNUSED(endJacobianPos, endStmtPos, endLLFByteDataPos, endLLFInfoDataPos);

        using Adjoint = AdjointVectorTraits::Gradient<AdjointVector>;

        typename Base::template VectorAccess<AdjointVector> vectorAccess(adjointVector);

        size_t curAdjointPos = startAdjointPos;

        while (curAdjointPos > endAdjointPos) CODI_Likely {
          curStmtPos -= 1;
          Config::ArgumentSize const argsSize = numberOfJacobians[curStmtPos];

          if (Config::StatementLowLevelFunctionTag == argsSize) CODI_Unlikely {
            Base::template callLowLevelFunction<LowLevelFunctionEntryCallKind::Reverse>(
                tape, false, curLLFByteDataPos, dataPtr, curLLFInfoDataPos, tokenPtr, dataSizePtr, &vectorAccess);
          } else if (Config::StatementInputTag == argsSize) CODI_Unlikely {
            // Do nothing.
          } else CODI_Likely {
            // No input value, perform regular statement evaluation.

            Adjoint const lhsAdjoint = adjointVector[curAdjointPos];  // We do not use the zero index, decrement of
                                                                      // curAdjointPos at the end of the loop.

            EventSystem<JacobianLinearTape>::notifyStatementEvaluateListeners(
                tape, (Identifier)curAdjointPos, GradientTraits::dim<Adjoint>(),
                GradientTraits::toArray(lhsAdjoint).data());

            if (Config::ReversalZeroesAdjoints) {
              adjointVector[curAdjointPos] = Adjoint();
            }

            Base::incrementAdjoints(adjointVector, lhsAdjoint, argsSize, curJacobianPos, rhsJacobians, rhsIdentifiers);
          }

          curAdjointPos -= 1;
        }
      }

      /// Passes the statement to the writer.
      template<typename TapeTypes>
      CODI_INLINE static void internalWriteTape(
          /* file interface pointer*/
          codi::TapeWriterInterface<TapeTypes>* writer,
          /* data from low level function byte data vector */
          size_t& curLLFByteDataPos, size_t const& endLLFByteDataPos, char* dataPtr,
          /* data from low level function info data vector */
          size_t& curLLFInfoDataPos, size_t const& endLLFInfoDataPos, Config::LowLevelFunctionToken* const tokenPtr,
          Config::LowLevelFunctionDataSize* const dataSizePtr,
          /* data from jacobianData */
          size_t& curJacobianPos, size_t const& endJacobianPos, typename Real::Real const* const rhsJacobians,
          typename Real::Identifier const* const rhsIdentifiers,
          /* data from statementData */
          size_t& curStmtPos, size_t const& endStmtPos, Config::ArgumentSize const* const numberOfJacobians,
          /* data from index handler */
          size_t const& startAdjointPos, size_t const& endAdjointPos) {
        CODI_UNUSED(endLLFByteDataPos, endLLFInfoDataPos, endJacobianPos, endStmtPos);

        ByteDataView dataView = {};
        LowLevelFunctionEntry<JacobianLinearTape, Real, Identifier> const* func = nullptr;

        typename Real::Identifier curLhsIdentifier;
        Config::ArgumentSize argsSize;
        size_t curAdjointPos = startAdjointPos;
        while (curAdjointPos < endAdjointPos) CODI_Likely {
          curAdjointPos += 1;
          curLhsIdentifier = curAdjointPos;
          argsSize = numberOfJacobians[curStmtPos];
          if (Config::StatementLowLevelFunctionTag == argsSize) CODI_Unlikely {
            Base::prepareLowLevelFunction(true, curLLFByteDataPos, dataPtr, curLLFInfoDataPos, tokenPtr, dataSizePtr,
                                          dataView, func);
            writer->writeLowLevelFunction(func, dataView);
          } else if (Config::StatementInputTag == argsSize) CODI_Unlikely {
            writer->writeStatement(curLhsIdentifier, curJacobianPos, rhsJacobians, rhsIdentifiers, argsSize);
          } else CODI_Likely {
            writer->writeStatement(curLhsIdentifier, curJacobianPos, rhsJacobians, rhsIdentifiers, argsSize);
            curJacobianPos += argsSize;
          }
          curStmtPos += 1;
        }
      }

    public:
      /*******************************************************************************/
      /// @name Functions from CustomIteratorTapeInterface
      /// @{

      using Base::iterateForward;
      /// \copydoc codi::CustomIteratorTapeInterface::iterateForward
      template<typename Callbacks>
      CODI_INLINE void iterateForward(Callbacks&& callbacks, Position start, Position end) {
        auto evalFunc =
            [&callbacks](
                /* data from call */
                JacobianLinearTape& tape,
                /* data from low level function byte data vector */
                size_t& curLLFByteDataPos, size_t const& endLLFByteDataPos, char* dataPtr,
                /* data from low level function info data vector */
                size_t& curLLFInfoDataPos, size_t const& endLLFInfoDataPos,
                Config::LowLevelFunctionToken* const tokenPtr, Config::LowLevelFunctionDataSize* const dataSizePtr,
                /* data from jacobian vector */
                size_t& curJacobianPos, size_t const& endJacobianPos, Real* const rhsJacobians,
                Identifier* const rhsIdentifiers,
                /* data from statement vector */
                size_t& curStmtPos, size_t const& endStmtPos, Config::ArgumentSize const* const numberOfJacobians,
                /* data from index handler */
                size_t const& startAdjointPos, size_t const& endAdjointPos) CODI_LAMBDA_INLINE {
              CODI_UNUSED(tape, endJacobianPos, endLLFByteDataPos, endLLFInfoDataPos, endStmtPos);

              ByteDataView dataView = {};
              LowLevelFunctionEntry<JacobianLinearTape, Real, Identifier> const* func = nullptr;

              size_t curAdjointPos = startAdjointPos;
              while (curAdjointPos < endAdjointPos) CODI_Likely {
                curAdjointPos += 1;

                Config::ArgumentSize argsSize = numberOfJacobians[curStmtPos];

                if (Config::StatementLowLevelFunctionTag == argsSize) CODI_Unlikely {
                  Base::prepareLowLevelFunction(true, curLLFByteDataPos, dataPtr, curLLFInfoDataPos, tokenPtr,
                                                dataSizePtr, dataView, func);
                  callbacks.handleLowLevelFunction(*func, dataView);
                } else CODI_Likely {
                  if (Config::StatementInputTag == argsSize) CODI_Unlikely {
                    argsSize = 0;
                  }

                  Identifier lhsIdentifier = curAdjointPos;
                  callbacks.handleStatement(lhsIdentifier, argsSize, &rhsJacobians[curJacobianPos],
                                            &rhsIdentifiers[curJacobianPos]);

                  codiAssert(lhsIdentifier ==
                             (Identifier)curAdjointPos);  // Lhs identifiers can not be edited in a linear tape.

                  curJacobianPos += argsSize;
                }

                curStmtPos += 1;
              }
            };

        Base::llfByteData.evaluateForward(start, end, evalFunc, *this);
      }

      using Base::iterateReverse;
      /// \copydoc codi::CustomIteratorTapeInterface::iterateReverse
      template<typename Callbacks>
      CODI_INLINE void iterateReverse(Callbacks&& callbacks, Position start, Position end) {
        auto evalFunc =
            [&callbacks](
                /* data from call */
                JacobianLinearTape& tape,
                /* data from low level function byte data vector */
                size_t& curLLFByteDataPos, size_t const& endLLFByteDataPos, char* dataPtr,
                /* data from low level function info data vector */
                size_t& curLLFInfoDataPos, size_t const& endLLFInfoDataPos,
                Config::LowLevelFunctionToken* const tokenPtr, Config::LowLevelFunctionDataSize* const dataSizePtr,
                /* data from jacobian vector */
                size_t& curJacobianPos, size_t const& endJacobianPos, Real* const rhsJacobians,
                Identifier* const rhsIdentifiers,
                /* data from statement vector */
                size_t& curStmtPos, size_t const& endStmtPos, Config::ArgumentSize const* const numberOfJacobians,
                /* data from index handler */
                size_t const& startAdjointPos, size_t const& endAdjointPos) CODI_LAMBDA_INLINE {
              CODI_UNUSED(tape, endJacobianPos, endLLFByteDataPos, endLLFInfoDataPos, endStmtPos);

              ByteDataView dataView = {};
              LowLevelFunctionEntry<JacobianLinearTape, Real, Identifier> const* func = nullptr;

              size_t curAdjointPos = startAdjointPos;
              while (curAdjointPos > endAdjointPos) CODI_Likely {
                curStmtPos -= 1;

                Config::ArgumentSize argsSize = numberOfJacobians[curStmtPos];

                if (Config::StatementLowLevelFunctionTag == argsSize) CODI_Unlikely {
                  Base::prepareLowLevelFunction(false, curLLFByteDataPos, dataPtr, curLLFInfoDataPos, tokenPtr,
                                                dataSizePtr, dataView, func);
                  callbacks.handleLowLevelFunction(*func, dataView);
                } else CODI_Likely {
                  if (Config::StatementInputTag == argsSize) CODI_Unlikely {
                    argsSize = 0;
                  }

                  curJacobianPos -= argsSize;

                  Identifier lhsIdentifier = curAdjointPos;
                  callbacks.handleStatement(lhsIdentifier, argsSize, &rhsJacobians[curJacobianPos],
                                            &rhsIdentifiers[curJacobianPos]);

                  codiAssert(lhsIdentifier ==
                             (Identifier)curAdjointPos);  // Lhs identifiers can not be edited in a linear tape.
                }

                curAdjointPos -= 1;
              }
            };

        Base::llfByteData.evaluateForward(start, end, evalFunc, *this);
      }

      /// @}
  };
}
