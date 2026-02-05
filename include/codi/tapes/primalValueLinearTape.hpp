/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2026 Chair for Scientific Computing (SciComp), RPTU University Kaiserslautern-Landau
 * Homepage: http://scicomp.rptu.de
 * Contact:  Prof. Nicolas R. Gauger (codi@scicomp.uni-kl.de)
 *
 * Lead developers: Max Sagebaum, Johannes Blühdorn (SciComp, RPTU University Kaiserslautern-Landau)
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
 *  - SciComp, RPTU University Kaiserslautern-Landau:
 *    - Max Sagebaum
 *    - Johannes Blühdorn
 *    - Former members:
 *      - Tim Albring
 */
#pragma once

#include <algorithm>
#include <functional>
#include <type_traits>

#include "../config.h"
#include "../expressions/lhsExpressionInterface.hpp"
#include "../expressions/logic/compileTimeTraversalLogic.hpp"
#include "../expressions/logic/constructStaticContext.hpp"
#include "../expressions/logic/traversalLogic.hpp"
#include "../misc/macros.hpp"
#include "../misc/memberStore.hpp"
#include "../traits/expressionTraits.hpp"
#include "data/chunk.hpp"
#include "indices/indexManagerInterface.hpp"
#include "primalValueBaseTape.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Final implementation for a primal value tape with a linear index management.
   *
   * This class implements the interface methods from the PrimalValueBaseTape.
   *
   * @tparam T_TapeTypes  JacobianTapeTypes definition.
   */
  template<typename T_TapeTypes>
  struct PrimalValueLinearTape : public PrimalValueBaseTape<T_TapeTypes, PrimalValueLinearTape<T_TapeTypes>> {
    public:

      using TapeTypes =
          CODI_DD(T_TapeTypes,
                  CODI_T(PrimalValueTapeTypes<double, double, IndexManagerInterface<int>, StatementEvaluatorInterface,
                                              DefaultChunkedData>));  ///< See PrimalValueLinearTape.

      using Base = PrimalValueBaseTape<T_TapeTypes, PrimalValueLinearTape<T_TapeTypes>>;  ///< Base class abbreviation.
      friend Base;  ///< Allow the base class to call protected and private methods.

      using Real = typename TapeTypes::Real;                              ///< See TapeTypesInterface.
      using Gradient = typename TapeTypes::Gradient;                      ///< See TapeTypesInterface.
      using IndexManager = typename TapeTypes::IndexManager;              ///< See PrimalValueTapeTypes.
      using Identifier = typename TapeTypes::Identifier;                  ///< See TapeTypesInterface.
      using PassiveReal = RealTraits::PassiveReal<Real>;                  ///< Basic computation type.
      using StatementEvaluator = typename TapeTypes::StatementEvaluator;  ///< See PrimalValueTapeTypes.
      using EvalHandle = typename TapeTypes::EvalHandle;                  ///< See PrimalValueTapeTypes.
      using Position = typename Base::Position;                           ///< See TapeTypesInterface.

      using StatementDataPointers = typename Base::StatementDataPointers;  ///< Defined in PrimalValueBaseTape.

      template<typename T>
      using StackArray = typename Base::template StackArray<T>;  ///< See PrimalValueBaseTape.

      /// Constructor
      PrimalValueLinearTape() : Base() {}

      using Base::clearAdjoints;

      /// \copydoc codi::PositionalEvaluationTapeInterface::clearAdjoints
      /// <br> Implementation: Automatic adjoints management has no effect. Primal value tapes do not implement adjoints
      /// locking.
      void clearAdjoints(Position const& start, Position const& end,
                         AdjointsManagement adjointsManagement = AdjointsManagement::Automatic) {
        CODI_UNUSED(adjointsManagement);

        using IndexPosition = CODI_DD(typename IndexManager::Position, int);
        IndexPosition startIndex = this->llfByteData.template extractPosition<IndexPosition>(start);
        IndexPosition endIndex = this->llfByteData.template extractPosition<IndexPosition>(end);

        startIndex = std::min(startIndex, (IndexPosition)this->adjoints.size() - 1);
        endIndex = std::min(endIndex, (IndexPosition)this->adjoints.size() - 1);

        for (IndexPosition curPos = endIndex + 1; curPos <= startIndex; curPos += 1) {
          this->adjoints[curPos] = Gradient();
        }
      }

      /// \copydoc codi::CustomAdjointVectorEvaluationTapeInterface::clearCustomAdjoints
      template<typename AdjointVector>
      void clearCustomAdjoints(Position const& start, Position const& end, AdjointVector data) {
        using IndexPosition = CODI_DD(typename IndexManager::Position, int);
        IndexPosition startIndex = this->llfByteData.template extractPosition<IndexPosition>(start);
        IndexPosition endIndex = this->llfByteData.template extractPosition<IndexPosition>(end);

        for (IndexPosition curPos = endIndex + 1; curPos <= startIndex; curPos += 1) {
          data[curPos] = AdjointVectorTraits::Gradient<AdjointVector>();
        }
      }

    protected:

      /// \copydoc codi::PrimalValueBaseTape::internalEvaluateForward_EvalStatements
      CODI_INLINE static void internalEvaluateForward_EvalStatements(
          /* data from call */
          PrimalValueLinearTape& tape, Real* primalVector, ADJOINT_VECTOR_TYPE* adjointVector,
          /* data from low level function byte data vector */
          size_t& curLLFByteDataPos, size_t const& endLLFByteDataPos, char* dataPtr,
          /* data from low level function info data vector */
          size_t& curLLFInfoDataPos, size_t const& endLLFInfoDataPos, Config::LowLevelFunctionToken* const tokenPtr,
          Config::LowLevelFunctionDataSize* const dataSizePtr,
          /* data from statementByteData */
          size_t& curStatementBytePos, size_t const& endStatementBytePos, char* stmtDataPtr,
          /* data from statementData */
          size_t& curStatementPos, size_t const& endStatementPos,
          Config::ArgumentSize const* const numberOfPassiveArguments, EvalHandle const* const stmtEvalHandle,
          Config::LowLevelFunctionDataSize* const stmtByteSize,
          /* data from index handler */
          size_t const& startAdjointPos, size_t const& endAdjointPos) {
        CODI_UNUSED(endLLFByteDataPos, endLLFInfoDataPos, endStatementBytePos, endAdjointPos);

        size_t curAdjointPos = startAdjointPos;
        StackArray<Real> lhsPrimals = {};
        StackArray<Gradient> lhsTangents = {};

#if !CODI_VariableAdjointInterfaceInPrimalTapes
        typename Base::template VectorAccess<Gradient*> vectorAccess(adjointVector, primalVector);
#endif

        while (curStatementPos < endStatementPos) CODI_Likely {
          Config::ArgumentSize nPassiveValues = numberOfPassiveArguments[curStatementPos];

          if (Config::StatementLowLevelFunctionTag == nPassiveValues) CODI_Unlikely {
            Base::template callLowLevelFunction<LowLevelFunctionEntryCallKind::Forward>(
                tape, true, curLLFByteDataPos, dataPtr, curLLFInfoDataPos, tokenPtr, dataSizePtr,
#if CODI_VariableAdjointInterfaceInPrimalTapes
                adjointVector
#else
                &vectorAccess
#endif
            );
          } else if (Config::StatementInputTag == nPassiveValues) CODI_Unlikely {
            curAdjointPos += 1;
          } else CODI_Likely {
            StatementEvaluator::template call<StatementCall::Forward, PrimalValueLinearTape>(
                stmtEvalHandle[curStatementPos], tape, lhsPrimals.data(), lhsTangents.data(), primalVector,
                adjointVector, curAdjointPos, nPassiveValues, &stmtDataPtr[curStatementBytePos]);

            curStatementBytePos += stmtByteSize[curStatementPos];
          }

          curStatementPos += 1;
        }
      }

      /// \copydoc codi::PrimalValueBaseTape::internalEvaluatePrimal_EvalStatements
      CODI_INLINE static void internalEvaluatePrimal_EvalStatements(
          /* data from call */
          PrimalValueLinearTape& tape, Real* primalVector,
          /* data from low level function byte data vector */
          size_t& curLLFByteDataPos, size_t const& endLLFByteDataPos, char* dataPtr,
          /* data from low level function info data vector */
          size_t& curLLFInfoDataPos, size_t const& endLLFInfoDataPos, Config::LowLevelFunctionToken* const tokenPtr,
          Config::LowLevelFunctionDataSize* const dataSizePtr,
          /* data from statementByteData */
          size_t& curStatementBytePos, size_t const& endStatementBytePos, char* stmtDataPtr,
          /* data from statementData */
          size_t& curStatementPos, size_t const& endStatementPos,
          Config::ArgumentSize const* const numberOfPassiveArguments, EvalHandle const* const stmtEvalHandle,
          Config::LowLevelFunctionDataSize* const stmtByteSize,
          /* data from index handler */
          size_t const& startAdjointPos, size_t const& endAdjointPos) {
        CODI_UNUSED(endLLFByteDataPos, endLLFInfoDataPos, endStatementBytePos, endAdjointPos);

        size_t curAdjointPos = startAdjointPos;
        StackArray<Real> lhsPrimals = {};

        typename Base::template VectorAccess<Gradient*> vectorAccess(nullptr, primalVector);

        while (curStatementPos < endStatementPos) CODI_Likely {
          Config::ArgumentSize nPassiveValues = numberOfPassiveArguments[curStatementPos];

          if (Config::StatementLowLevelFunctionTag == nPassiveValues) CODI_Unlikely {
            Base::template callLowLevelFunction<LowLevelFunctionEntryCallKind::Primal>(
                tape, true, curLLFByteDataPos, dataPtr, curLLFInfoDataPos, tokenPtr, dataSizePtr, &vectorAccess);
          } else if (Config::StatementInputTag == nPassiveValues) CODI_Unlikely {
            curAdjointPos += 1;
          } else CODI_Likely {
            StatementEvaluator::template call<StatementCall::Primal, PrimalValueLinearTape>(
                stmtEvalHandle[curStatementPos], tape, lhsPrimals.data(), primalVector, curAdjointPos, nPassiveValues,
                &stmtDataPtr[curStatementBytePos]);

            curStatementBytePos += stmtByteSize[curStatementPos];
          }

          curStatementPos += 1;
        }
      }

      /// \copydoc codi::PrimalValueBaseTape::internalEvaluateReverse_EvalStatements
      CODI_INLINE static void internalEvaluateReverse_EvalStatements(
          /* data from call */
          PrimalValueLinearTape& tape, Real* primalVector, ADJOINT_VECTOR_TYPE* adjointVector,
          /* data from low level function byte data vector */
          size_t& curLLFByteDataPos, size_t const& endLLFByteDataPos, char* dataPtr,
          /* data from low level function info data vector */
          size_t& curLLFInfoDataPos, size_t const& endLLFInfoDataPos, Config::LowLevelFunctionToken* const tokenPtr,
          Config::LowLevelFunctionDataSize* const dataSizePtr,
          /* data from statementByteData */
          size_t& curStatementBytePos, size_t const& endStatementBytePos, char* stmtDataPtr,
          /* data from statementData */
          size_t& curStatementPos, size_t const& endStatementPos,
          Config::ArgumentSize const* const numberOfPassiveArguments, EvalHandle const* const stmtEvalHandle,
          Config::LowLevelFunctionDataSize* const stmtByteSize,
          /* data from index handler */
          size_t const& startAdjointPos, size_t const& endAdjointPos) {
        CODI_UNUSED(endLLFByteDataPos, endLLFInfoDataPos, endStatementBytePos, endAdjointPos);

        size_t curAdjointPos = startAdjointPos;
        StackArray<Gradient> lhsAdjoints = {};

#if !CODI_VariableAdjointInterfaceInPrimalTapes
        typename Base::template VectorAccess<Gradient*> vectorAccess(adjointVector, primalVector);
#endif

        while (curStatementPos > endStatementPos) CODI_Likely {
          curStatementPos -= 1;

          Config::ArgumentSize nPassiveValues = numberOfPassiveArguments[curStatementPos];

          if (Config::StatementLowLevelFunctionTag == nPassiveValues) CODI_Unlikely {
            Base::template callLowLevelFunction<LowLevelFunctionEntryCallKind::Reverse>(
                tape, false, curLLFByteDataPos, dataPtr, curLLFInfoDataPos, tokenPtr, dataSizePtr,
#if CODI_VariableAdjointInterfaceInPrimalTapes
                adjointVector
#else
                &vectorAccess
#endif
            );
          } else if (Config::StatementInputTag == nPassiveValues) CODI_Unlikely {
            curAdjointPos -= 1;
          } else CODI_Likely {
            curStatementBytePos -= stmtByteSize[curStatementPos];

            StatementEvaluator::template call<StatementCall::Reverse, PrimalValueLinearTape>(
                stmtEvalHandle[curStatementPos], tape, lhsAdjoints.data(), primalVector, adjointVector, curAdjointPos,
                nPassiveValues, &stmtDataPtr[curStatementBytePos]);
          }
        }
      }

      /// Passes the statement information and the stmtEvalHandle to the writer.
      template<typename TapeTypes>
      CODI_INLINE static void internalWriteTape(
          /* data from call */
          Real* primalVector,
          /* file interface pointer*/
          codi::TapeWriterInterface<TapeTypes>* writer,
          /* data from low level function byte data vector */
          size_t& curLLFByteDataPos, size_t const& endLLFByteDataPos, char* dataPtr,
          /* data from low level function info data vector */
          size_t& curLLFInfoDataPos, size_t const& endLLFInfoDataPos, Config::LowLevelFunctionToken* const tokenPtr,
          Config::LowLevelFunctionDataSize* const dataSizePtr,
          /* data from statementByteData */
          size_t& curStatementBytePos, size_t const& endStatementBytePos, char* stmtDataPtr,
          /* data from statementData */
          size_t& curStatementPos, size_t const& endStatementPos,
          Config::ArgumentSize const* const numberOfPassiveArguments, EvalHandle const* const stmtEvalHandle,
          Config::LowLevelFunctionDataSize* const stmtByteSize,
          /* data from index handler */
          size_t const& startAdjointPos, size_t const& endAdjointPos) {
        CODI_UNUSED(endLLFByteDataPos, endLLFInfoDataPos, endStatementBytePos, endAdjointPos);

        size_t curAdjointPos = startAdjointPos;
        StackArray<Identifier> lhsIdentifiers;

        ByteDataView dataView = {};
        LowLevelFunctionEntry<PrimalValueLinearTape, Real, Identifier> const* func = nullptr;

        while (curStatementPos < endStatementPos) {
          Config::ArgumentSize nPassiveValues = numberOfPassiveArguments[curStatementPos];

          if (Config::StatementLowLevelFunctionTag == nPassiveValues) CODI_Unlikely {
            Base::prepareLowLevelFunction(true, curLLFByteDataPos, dataPtr, curLLFInfoDataPos, tokenPtr, dataSizePtr,
                                          dataView, func);
            writer->writeLowLevelFunction(func, dataView);
          } else CODI_Likely {
            WriteInfo writeInfo;
            StatementEvaluator::template call<StatementCall::WriteInformation, PrimalValueLinearTape>(
                stmtEvalHandle[curStatementPos], writeInfo, primalVector, nPassiveValues,
                &stmtDataPtr[curStatementBytePos]);

            StatementDataPointers pointers = {};
            pointers.populate(writeInfo.numberOfOutputArguments, writeInfo.numberOfActiveArguments, nPassiveValues,
                              writeInfo.numberOfConstantArguments, &stmtDataPtr[curStatementBytePos]);

            Real const* lhsPrimalValues = &primalVector[curAdjointPos + 1];
            for (size_t i = 0; i < writeInfo.numberOfOutputArguments; i += 1) {
              curAdjointPos += 1;
              lhsIdentifiers[i] = curAdjointPos;
            }

            writer->writeStatement(writeInfo, lhsIdentifiers.data(), lhsPrimalValues, nPassiveValues,
                                   pointers.rhsIdentifiers, pointers.passiveValues, pointers.constantValues,
                                   stmtEvalHandle[curStatementPos]);

            curStatementBytePos += stmtByteSize[curStatementPos];
          }

          curStatementPos += 1;
        }
      }

      /// \copydoc codi::PrimalValueBaseTape::internalResetPrimalValues
      /// Empty implementation; primal values are not overwritten with linear index management.
      CODI_INLINE void internalResetPrimalValues(Position const& pos) {
        CODI_UNUSED(pos);

        // Nothing to do.
      }

    public:
      /// \copydoc codi::PrimalValueBaseTape::revertPrimals
      /// Empty implementation; primal values are not overwritten with linear index management.
      void revertPrimals(Position const& pos) {
        CODI_UNUSED(pos);

        // Primal values do not need to be reset.
      }

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
                PrimalValueLinearTape& tape,
                /* data from low level function byte data vector */
                size_t& curLLFByteDataPos, size_t const& endLLFByteDataPos, char* dataPtr,
                /* data from low level function info data vector */
                size_t& curLLFInfoDataPos, size_t const& endLLFInfoDataPos,
                Config::LowLevelFunctionToken* const tokenPtr, Config::LowLevelFunctionDataSize* const dataSizePtr,
                /* data from statementByteData */
                size_t& curStatementBytePos, size_t const& endStatementBytePos, char* stmtDataPtr,
                /* data from statementData */
                size_t& curStatementPos, size_t const& endStatementPos,
                Config::ArgumentSize const* const numberOfPassiveArguments, EvalHandle const* const stmtEvalHandle,
                Config::LowLevelFunctionDataSize* const stmtByteSize,
                /* data from index handler */
                size_t const& startAdjointPos, size_t const& endAdjointPos) CODI_LAMBDA_INLINE {
              CODI_UNUSED(tape, endLLFByteDataPos, endLLFInfoDataPos, endStatementBytePos, endAdjointPos);

              size_t curAdjointPos = startAdjointPos;
              ByteDataView dataView = {};
              LowLevelFunctionEntry<PrimalValueLinearTape, Real, Identifier> const* func = nullptr;

              while (curStatementPos < endStatementPos) CODI_Likely {
                Config::ArgumentSize nPassiveValues = numberOfPassiveArguments[curStatementPos];

                if (Config::StatementLowLevelFunctionTag == nPassiveValues) CODI_Unlikely {
                  Base::prepareLowLevelFunction(true, curLLFByteDataPos, dataPtr, curLLFInfoDataPos, tokenPtr,
                                                dataSizePtr, dataView, func);
                  callbacks.handleLowLevelFunction(*func, dataView);
                } else CODI_Likely {
                  if (Config::StatementInputTag == nPassiveValues) CODI_Unlikely {
                    nPassiveValues = 0;
                  }
                  callbacks.handleStatement(stmtEvalHandle[curStatementPos], nPassiveValues, curAdjointPos,
                                            &stmtDataPtr[curStatementBytePos]);

                  curStatementBytePos += stmtByteSize[curStatementPos];
                }

                curStatementPos += 1;
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
                PrimalValueLinearTape& tape,
                /* data from low level function byte data vector */
                size_t& curLLFByteDataPos, size_t const& endLLFByteDataPos, char* dataPtr,
                /* data from low level function info data vector */
                size_t& curLLFInfoDataPos, size_t const& endLLFInfoDataPos,
                Config::LowLevelFunctionToken* const tokenPtr, Config::LowLevelFunctionDataSize* const dataSizePtr,
                /* data from statementByteData */
                size_t& curStatementBytePos, size_t const& endStatementBytePos, char* stmtDataPtr,
                /* data from statementData */
                size_t& curStatementPos, size_t const& endStatementPos,
                Config::ArgumentSize const* const numberOfPassiveArguments, EvalHandle const* const stmtEvalHandle,
                Config::LowLevelFunctionDataSize* const stmtByteSize,
                /* data from index handler */
                size_t const& startAdjointPos, size_t const& endAdjointPos) CODI_LAMBDA_INLINE {
              CODI_UNUSED(tape, endLLFByteDataPos, endLLFInfoDataPos, endStatementBytePos, endAdjointPos);

              size_t curAdjointPos = startAdjointPos;
              ByteDataView dataView = {};
              LowLevelFunctionEntry<PrimalValueLinearTape, Real, Identifier> const* func = nullptr;

              while (curStatementPos > endStatementPos) CODI_Likely {
                curStatementPos -= 1;

                Config::ArgumentSize nPassiveValues = numberOfPassiveArguments[curStatementPos];

                if (Config::StatementLowLevelFunctionTag == nPassiveValues) CODI_Unlikely {
                  Base::prepareLowLevelFunction(false, curLLFByteDataPos, dataPtr, curLLFInfoDataPos, tokenPtr,
                                                dataSizePtr, dataView, func);
                  callbacks.handleLowLevelFunction(*func, dataView);
                } else CODI_Likely {
                  curStatementBytePos -= stmtByteSize[curStatementPos];

                  if (Config::StatementInputTag == nPassiveValues) CODI_Unlikely {
                    nPassiveValues = 0;
                  }
                  callbacks.handleStatement(stmtEvalHandle[curStatementPos], nPassiveValues, curAdjointPos,
                                            &stmtDataPtr[curStatementBytePos]);
                }
              }
            };

        Base::llfByteData.evaluateReverse(start, end, evalFunc, *this);
      }
  };
}
