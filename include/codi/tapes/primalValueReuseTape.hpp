/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2026 Chair for Scientific Computing (SciComp), University of Kaiserslautern-Landau
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
   * @brief Final implementation for a primal value tape with a reuse index management.
   *
   * This class implements the interface methods from the PrimalValueBaseTape.
   *
   * @tparam T_TapeTypes  JacobianTapeTypes definition.
   */
  template<typename T_TapeTypes>
  struct PrimalValueReuseTape : public PrimalValueBaseTape<T_TapeTypes, PrimalValueReuseTape<T_TapeTypes>> {
    public:

      using TapeTypes =
          CODI_DD(T_TapeTypes,
                  CODI_T(PrimalValueTapeTypes<double, double, IndexManagerInterface<int>, StatementEvaluatorInterface,
                                              DefaultChunkedData>));  ///< See PrimalValueReuseTape.

      using Base = PrimalValueBaseTape<T_TapeTypes, PrimalValueReuseTape<T_TapeTypes>>;  ///< Base class abbreviation.
      friend Base;  ///< Allow the base class to call protected and private methods.

      using Real = typename TapeTypes::Real;                              ///< See TapeTypesInterface.
      using Gradient = typename TapeTypes::Gradient;                      ///< See TapeTypesInterface.
      using Identifier = typename TapeTypes::Identifier;                  ///< See TapeTypesInterface.
      using PassiveReal = RealTraits::PassiveReal<Real>;                  ///< Basic computation type.
      using StatementEvaluator = typename TapeTypes::StatementEvaluator;  ///< See PrimalValueTapeTypes.
      using EvalHandle = typename TapeTypes::EvalHandle;                  ///< See PrimalValueTapeTypes.
      using Position = typename Base::Position;                           ///< See TapeTypesInterface.

      using StatementData = typename TapeTypes::StatementData;  ///< See PrimalValueTapeTypes.

      using StatementDataPointers = typename Base::StatementDataPointers;  ///< Defined in PrimalValueBaseTape.

      template<typename T>
      using StackArray = typename Base::template StackArray<T>;  ///< See PrimalValueBaseTape.

      /// Constructor
      PrimalValueReuseTape() : Base() {}

      using Base::clearAdjoints;

    private:

      CODI_INLINE static void internalClearAdjoints_EvalStatements(
          /* data from call */
          ADJOINT_VECTOR_TYPE* adjointVector,
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
          Config::LowLevelFunctionDataSize* const stmtByteSize) {
        CODI_UNUSED(endLLFByteDataPos, endLLFInfoDataPos, endStatementBytePos);

        while (curStatementPos < endStatementPos) CODI_Likely {
          Config::ArgumentSize nPassiveValues = numberOfPassiveArguments[curStatementPos];

          if (Config::StatementLowLevelFunctionTag == nPassiveValues) CODI_Unlikely {
            Base::skipLowLevelFunction(true, curLLFByteDataPos, dataPtr, curLLFInfoDataPos, tokenPtr, dataSizePtr);
          } else CODI_Likely {
            StatementEvaluator::template call<StatementCall::ClearAdjoints, PrimalValueReuseTape>(
                stmtEvalHandle[curStatementPos], adjointVector, nPassiveValues, &stmtDataPtr[curStatementBytePos]);

            curStatementBytePos += stmtByteSize[curStatementPos];
          }

          curStatementPos += 1;
        }
      }

    public:

      /// \copydoc codi::PositionalEvaluationTapeInterface::clearAdjoints
      /// <br> Implementation: Automatic adjoints management has no effect. Primal value tapes do not implement adjoints
      /// locking.
      void clearAdjoints(Position const& start, Position const& end,
                         AdjointsManagement adjointsManagement = AdjointsManagement::Automatic) {
        CODI_UNUSED(adjointsManagement);

        clearCustomAdjoints(start, end, this->adjoints.data());
      }

      /// \copydoc codi::CustomAdjointVectorEvaluationTapeInterface::clearCustomAdjoints
      template<typename AdjointVector>
      void clearCustomAdjoints(Position const& start, Position const& end, AdjointVector data) {
        CODI_STATIC_ASSERT(
            Config::VariableAdjointInterfaceInPrimalTapes ||
                CODI_T(std::is_same<typename std::remove_reference<AdjointVector>::type, Gradient*>::value),
            "Please enable 'CODI_VariableAdjointInterfaceInPrimalTapes' in order"
            " to use custom adjoint vectors in the primal value tapes.");

        typename Base::template VectorAccess<AdjointVector> vectorAccess(data, this->primals.data());

        ADJOINT_VECTOR_TYPE* dataVector = Base::template selectAdjointVector<AdjointVector>(&vectorAccess, data);

        this->llfByteData.evaluateForward(end, start, internalClearAdjoints_EvalStatements, dataVector);
      }

    protected:

      /// \copydoc codi::PrimalValueBaseTape::internalEvaluateForward_EvalStatements
      CODI_INLINE static void internalEvaluateForward_EvalStatements(
          /* data from call */
          PrimalValueReuseTape& tape, Real* primalVector, ADJOINT_VECTOR_TYPE* adjointVector,
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
          Config::LowLevelFunctionDataSize* const stmtByteSize) {
        CODI_UNUSED(endLLFByteDataPos, endLLFInfoDataPos, endStatementBytePos);

#if !CODI_VariableAdjointInterfaceInPrimalTapes
        typename Base::template VectorAccess<Gradient*> vectorAccess(adjointVector, primalVector);
#endif

        size_t linearAdjointPos = 0;  // Not accessed by the implementation, just a temporary.
        StackArray<Real> lhsPrimals = {};
        StackArray<Gradient> lhsTangents = {};

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
          } else CODI_Likely {
            StatementEvaluator::template call<StatementCall::Forward, PrimalValueReuseTape>(
                stmtEvalHandle[curStatementPos], tape, lhsPrimals.data(), lhsTangents.data(), primalVector,
                adjointVector, linearAdjointPos, nPassiveValues, &stmtDataPtr[curStatementBytePos]);

            curStatementBytePos += stmtByteSize[curStatementPos];
          }

          curStatementPos += 1;
        }
      }

      /// \copydoc codi::PrimalValueBaseTape::internalEvaluatePrimal_EvalStatements
      CODI_INLINE static void internalEvaluatePrimal_EvalStatements(
          /* data from call */
          PrimalValueReuseTape& tape, Real* primalVector,
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
          Config::LowLevelFunctionDataSize* const stmtByteSize) {
        CODI_UNUSED(endLLFByteDataPos, endLLFInfoDataPos, endStatementBytePos);

        typename Base::template VectorAccess<Gradient*> vectorAccess(nullptr, primalVector);
        StackArray<Real> lhsPrimals = {};

        size_t linearAdjointPos = 0;  // Not accessd by the implementation, just a temporary.

        while (curStatementPos < endStatementPos) CODI_Likely {
          Config::ArgumentSize nPassiveValues = numberOfPassiveArguments[curStatementPos];

          if (Config::StatementLowLevelFunctionTag == nPassiveValues) CODI_Unlikely {
            Base::template callLowLevelFunction<LowLevelFunctionEntryCallKind::Primal>(
                tape, true, curLLFByteDataPos, dataPtr, curLLFInfoDataPos, tokenPtr, dataSizePtr, &vectorAccess);
          } else CODI_Likely {
            StatementEvaluator::template call<StatementCall::Primal, PrimalValueReuseTape>(
                stmtEvalHandle[curStatementPos], tape, lhsPrimals.data(), primalVector, linearAdjointPos,
                numberOfPassiveArguments[curStatementPos], &stmtDataPtr[curStatementBytePos]);

            curStatementBytePos += stmtByteSize[curStatementPos];
          }

          curStatementPos += 1;
        }
      }

      /// \copydoc codi::PrimalValueBaseTape::internalEvaluateReverse_EvalStatements
      CODI_INLINE static void internalEvaluateReverse_EvalStatements(
          /* data from call */
          PrimalValueReuseTape& tape, Real* primalVector, ADJOINT_VECTOR_TYPE* adjointVector,
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
          Config::LowLevelFunctionDataSize* const stmtByteSize) {
        CODI_UNUSED(endLLFByteDataPos, endLLFInfoDataPos, endStatementBytePos);

#if !CODI_VariableAdjointInterfaceInPrimalTapes
        typename Base::template VectorAccess<Gradient*> vectorAccess(adjointVector, primalVector);
#endif

        size_t linearAdjointPos = 0;  // Not accessd by the implementation, just a temporary.
        StackArray<Gradient> lhsAdjoints = {};

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
          } else CODI_Likely {
            curStatementBytePos -= stmtByteSize[curStatementPos];

            StatementEvaluator::template call<StatementCall::Reverse, PrimalValueReuseTape>(
                stmtEvalHandle[curStatementPos], tape, lhsAdjoints.data(), primalVector, adjointVector,
                linearAdjointPos, numberOfPassiveArguments[curStatementPos], &stmtDataPtr[curStatementBytePos]);
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
          Config::LowLevelFunctionDataSize* const stmtByteSize) {
        CODI_UNUSED(endLLFByteDataPos, endLLFInfoDataPos, endStatementBytePos);

        ByteDataView dataView = {};
        LowLevelFunctionEntry<PrimalValueReuseTape, Real, Identifier> const* func = nullptr;

        while (curStatementPos < endStatementPos) {
          Config::ArgumentSize nPassiveValues = numberOfPassiveArguments[curStatementPos];

          if (Config::StatementLowLevelFunctionTag == nPassiveValues) CODI_Unlikely {
            Base::prepareLowLevelFunction(true, curLLFByteDataPos, dataPtr, curLLFInfoDataPos, tokenPtr, dataSizePtr,
                                          dataView, func);
            writer->writeLowLevelFunction(func, dataView);
          } else CODI_Likely {
            WriteInfo writeInfo;
            StatementEvaluator::template call<StatementCall::WriteInformation, PrimalValueReuseTape>(
                stmtEvalHandle[curStatementPos], writeInfo, primalVector, nPassiveValues,
                &stmtDataPtr[curStatementBytePos]);

            StatementDataPointers pointers = {};
            pointers.populate(writeInfo.numberOfOutputArguments, writeInfo.numberOfActiveArguments, nPassiveValues,
                              writeInfo.numberOfConstantArguments, &stmtDataPtr[curStatementBytePos]);

            writer->writeStatement(writeInfo, pointers.lhsIdentifiers, pointers.oldLhsValues, nPassiveValues,
                                   pointers.rhsIdentifiers, pointers.passiveValues, pointers.constantValues,
                                   stmtEvalHandle[curStatementPos]);

            curStatementBytePos += stmtByteSize[curStatementPos];
          }

          curStatementPos += 1;
        }
      }

    private:

      CODI_INLINE static void internalResetPrimals_EvalStatements(
          /* data from call */
          Real* primalVector,
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
          Config::LowLevelFunctionDataSize* const stmtByteSize) {
        CODI_UNUSED(endLLFByteDataPos, endLLFInfoDataPos, endStatementBytePos);

        while (curStatementPos > endStatementPos) CODI_Likely {
          curStatementPos -= 1;

          Config::ArgumentSize nPassiveValues = numberOfPassiveArguments[curStatementPos];

          if (Config::StatementLowLevelFunctionTag == nPassiveValues) CODI_Unlikely {
            Base::skipLowLevelFunction(false, curLLFByteDataPos, dataPtr, curLLFInfoDataPos, tokenPtr, dataSizePtr);
          } else CODI_Likely {
            curStatementBytePos -= stmtByteSize[curStatementPos];

            StatementEvaluator::template call<StatementCall::ResetPrimals, PrimalValueReuseTape>(
                stmtEvalHandle[curStatementPos], primalVector, nPassiveValues, &stmtDataPtr[curStatementBytePos]);
          }
        }
      }

    public:

      /// \copydoc codi::PrimalValueBaseTape::internalResetPrimalValues
      CODI_INLINE void internalResetPrimalValues(Position const& pos) {
        this->llfByteData.evaluateReverse(this->getPosition(), pos, internalResetPrimals_EvalStatements,
                                          this->primals.data());
      }

    public:
      /// \copydoc codi::PrimalEvaluationTapeInterface::revertPrimals
      void revertPrimals(Position const& pos) {
        internalResetPrimalValues(pos);
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
                PrimalValueReuseTape& tape,
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
                Config::LowLevelFunctionDataSize* const stmtByteSize) CODI_LAMBDA_INLINE {
              CODI_UNUSED(tape, endLLFByteDataPos, endLLFInfoDataPos, endStatementBytePos);

              size_t linearAdjointPos = 0;  // Not accessed by the implementation, just a temporary.
              ByteDataView dataView = {};
              LowLevelFunctionEntry<PrimalValueReuseTape, Real, Identifier> const* func = nullptr;

              while (curStatementPos < endStatementPos) CODI_Likely {
                Config::ArgumentSize nPassiveValues = numberOfPassiveArguments[curStatementPos];

                if (Config::StatementLowLevelFunctionTag == nPassiveValues) CODI_Unlikely {
                  Base::prepareLowLevelFunction(true, curLLFByteDataPos, dataPtr, curLLFInfoDataPos, tokenPtr,
                                                dataSizePtr, dataView, func);
                  callbacks.handleLowLevelFunction(*func, dataView);
                } else CODI_Likely {
                  callbacks.handleStatement(stmtEvalHandle[curStatementPos], nPassiveValues, linearAdjointPos,
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
                PrimalValueReuseTape& tape,
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
                Config::LowLevelFunctionDataSize* const stmtByteSize) CODI_LAMBDA_INLINE {
              CODI_UNUSED(tape, endLLFByteDataPos, endLLFInfoDataPos, endStatementBytePos);

              size_t linearAdjointPos = 0;  // Not accessed by the implementation, just a temporary.
              ByteDataView dataView = {};
              LowLevelFunctionEntry<PrimalValueReuseTape, Real, Identifier> const* func = nullptr;

              while (curStatementPos < endStatementPos) CODI_Likely {
                curStatementPos -= 1;
                Config::ArgumentSize nPassiveValues = numberOfPassiveArguments[curStatementPos];

                if (Config::StatementLowLevelFunctionTag == nPassiveValues) CODI_Unlikely {
                  Base::prepareLowLevelFunction(false, curLLFByteDataPos, dataPtr, curLLFInfoDataPos, tokenPtr,
                                                dataSizePtr, dataView, func);
                  callbacks.handleLowLevelFunction(*func, dataView);
                } else CODI_Likely {
                  curStatementBytePos -= stmtByteSize[curStatementPos];

                  callbacks.handleStatement(stmtEvalHandle[curStatementPos], nPassiveValues, linearAdjointPos,
                                            &stmtDataPtr[curStatementBytePos]);
                }
              }
            };

        Base::llfByteData.evaluateReverse(start, end, evalFunc, *this);
      }
  };
}
