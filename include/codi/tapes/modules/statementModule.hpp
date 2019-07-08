/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2019 Chair for Scientific Computing (SciComp), TU Kaiserslautern
 * Homepage: http://www.scicomp.uni-kl.de
 * Contact:  Prof. Nicolas R. Gauger (codi@scicomp.uni-kl.de)
 *
 * Lead developers: Max Sagebaum, Tim Albring (SciComp, TU Kaiserslautern)
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
 * Authors: Max Sagebaum, Tim Albring, (SciComp, TU Kaiserslautern)
 */

#pragma once

#include "../reverseTapeInterface.hpp"
#include "../../configure.h"
#include "../../tapeTypes.hpp"
#include "../../tools/jacobianSorter.hpp"
#include "../../tools/tapeValues.hpp"
#include "../../typeFunctions.hpp"

/**
 * @brief Global namespace for CoDiPack - Code Differentiation Package
 */
namespace codi {

  /**
   * The module defines the structures stmtVector.
   * The module defines the types StmtChildVector, StmtChildPosition, StmtVector, StmtChunk,
   * StmtPosition.
   *
   * It defines the methods store(Expr), store(const), store(User), printStmtStatistics from the TapeInterface and ReverseTapeInterface.
   *
   * It defines the methods setStatementChunkSize, getUsedStatementSize, evaluatePrimalStub, resizeStmt as interface functions for the
   * including class.
   *
   * @tparam TapeTypes  All the types for the tape. Including the calculation type and the vector types.
   * @tparam      Tape  The full tape implementation
   */
  template<typename TapeTypes, typename Tape>
  struct StatementModule : public virtual ReverseTapeInterface<typename TapeTypes::Real, typename TapeTypes::Index, typename TapeTypes::GradientValue, Tape, typename TapeTypes::Position > {

    private:

    // ----------------------------------------------------------------------
    // All definitions of the module
    // ----------------------------------------------------------------------

      CODI_INLINE_REVERSE_TAPE_TYPES(TapeTypes::BaseTypes)

      /** @brief The chunk vector for the statement data. */
      typedef typename TapeTypes::StatementVector StmtVector;

      /** @brief The child vector for the statement data vector. */
      typedef typename StmtVector::NestedVectorType StmtChildVector;

      /** @brief The position type of the statement child vector */
      typedef typename StmtChildVector::Position StmtChildPosition;

      /** @brief The data for each statement. */
      typedef typename StmtVector::ChunkType StmtChunk;

      /** @brief The position type of the statement module. */
      typedef typename StmtVector::Position StmtPosition;

      /** @brief The global position for the tape */
      typedef typename TapeTypes::Position Position;

      /** @brief Forward the GradientData from the tape types. */
      typedef typename TapeTypes::GradientData GradientData;

      /**
       * @brief Cast this class to the full.
       *
       * The full type is able to access all functions from the tape and other modules.
       *
       * @return  The full tape implemenation.
       */
      Tape& cast() {
        return *static_cast<Tape*>(this);
      }

    protected:

#if CODI_EnableCombineJacobianArguments
      /** @brief Helper structure to buffer the argument data and add arguments with the same identifier */
      JacobianSorter<Real, GradientData> insertData;
#endif

      /** @brief The data for the statements. */
      StmtVector stmtVector;

    private:

    public:

      StatementModule() :
#if CODI_EnableCombineJacobianArguments
        insertData(),
#endif
        stmtVector(DefaultChunkSize)
      {}

    protected:

      /**
       * @brief Adds statistics about the statements.
       *
       * Adds the number of chunks, the total number of statements, the
       * allocated memory and the used memory.
       *
       * @param[in,out] values  The values where the information is added to.
       */
      void addStmtValues(TapeValues& values) const {
        values.addSection("Statements");
        values.addStreamData(stmtVector);
      }

      /**
       * @brief Initialize the StatementModule.
       *
       * Called after all members of the tape have been initialized.
       *
       * @param[in,out] childVector  The child vector for the statement vector.
       */
      void initStmtModule(StmtChildVector* childVector) {
        stmtVector.setNested(childVector);
      }

      /**
       * @brief Perform all actions to compute the Jacobians of the expression and store them in the jacobian vector.
       *
       * @param[in] rhs  The right hand side expression
       *
       * @tparam Rhs  The expression of the right hand side.
       */
      template<typename Rhs>
      CODI_INLINE size_t addJacobianEntries(const Rhs& rhs) {

        size_t startSize = cast().jacobiVector.getChunkPosition();

#if !CODI_EnableCombineJacobianArguments
        // If enabled insertData is defined as a member
        auto& insertData = cast().jacobiVector;
#endif

        // Push the regular Jacobian arguments
        rhs.template calcGradient(insertData);

        // Push Jacobians from ReferencReal arguments
        rhs.template pushLazyJacobies(insertData);

        // Store the Jacobians if the cobine optimization is enabled
#if CODI_EnableCombineJacobianArguments
        insertData.storeData(cast().jacobiVector);
#endif
        return cast().jacobiVector.getChunkPosition() - startSize;
      }

    // ----------------------------------------------------------------------
    // Protected function for the communication with the including class
    // ----------------------------------------------------------------------

      /**
       * @brief Resize the statement data.
       *
       * Ensure that enough size is allocated such that dataSize number of items
       * can be stored.
       *
       * @param[in] statementSize  The size that should be allocated for the statement data.
       */
      void resizeStmt(const size_t& statementSize) {
        stmtVector.resize(statementSize);
      }

      /**
       * @brief Stub for no primal evaluation, since it is not supported for Jacobian tapes.
       *
       * @param[in] start  The starting position for the forward evaluation.
       * @param[in]   end  The ending position for the forward evaluation.
       */
      CODI_INLINE void evaluatePrimalInternal(const Position& start, const Position& end) {
        CODI_UNUSED(start);
        CODI_UNUSED(end);
      }

    public:

    // ----------------------------------------------------------------------
    // Public function from the TapeInterface and ReverseTapeInterface
    // ----------------------------------------------------------------------

      /**
       * @brief Set the size of the statement data chunks.
       *
       * @param[in] statementChunkSize The new size for the statement data chunks.
       */
      void setStatementChunkSize(const size_t& statementChunkSize) {
        stmtVector.setChunkSize(statementChunkSize);
      }


      /**
       * @brief Store the jacobies of the statement on the tape.
       *
       * The jacobies and the indices of the rhs expression are
       * stored on the tape. Also the number of active variables
       * is stored in the statement vector.
       *
       * The gradient data of the lhs will get a new index.
       * The primal value of the lhs is set to the primal value of the rhs.
       *
       * @param[out]   lhsValue    The primal value of the lhs. This value is set to the value
       *                           of the right hand side.
       * @param[out]   lhsIndex    The gradient data of the lhs. The index will be updated.
       * @param[in]         rhs    The right hand side expression of the assignment.
       *
       * @tparam Rhs The expression on the rhs of the statement.
       */
      template<typename Rhs>
      CODI_INLINE void store(Real& lhsValue, Index& lhsIndex, const Rhs& rhs) {
        Tape& tape = cast();

        static_assert(ExpressionTraits<Rhs>::maxActiveVariables < MaxStatementIntSize, "Expression with to many arguments.");

        ENABLE_CHECK (OptTapeActivity, tape.active){
          stmtVector.reserveItems(1);
          tape.jacobiVector.reserveItems(ExpressionTraits<Rhs>::maxActiveVariables);
          /* first store the size of the current stack position and evaluate the
           rhs expression. If there was an active variable on the rhs, update
           the index of the lhs */
          size_t activeVariables = addJacobianEntries(rhs);
          ENABLE_CHECK(OptCheckEmptyStatements, 0 != activeVariables) {

            tape.indexHandler.assignIndex(lhsIndex);
            tape.pushStmtData((StatementInt)activeVariables, lhsIndex);

#if CODI_AdjointHandle_Jacobi
            Real* jacobies = NULL;
            Index* rhsIndices = NULL;

            JACOBI_VECTOR_NAME.getDataPointer(jacobies, rhsIndices);
            jacobies -= activeVariables;
            rhsIndices -= activeVariables;

            handleAdjointOperation(rhs.getValue(), lhsIndex, jacobies, rhsIndices, activeVariables);
#endif
          } else {
            tape.indexHandler.freeIndex(lhsIndex);
          }
        } else {
          tape.indexHandler.freeIndex(lhsIndex);
        }
        /* now set the value of the lhs */
        lhsValue = rhs.getValue();
      }

      /**
       * @brief Optimization for a passive value on the rhs. The lhs index is set to zero.
       *
       * No data is stored in this method.
       *
       * The primal value of the lhs is set to the primal value of the rhs.
       *
       * @param[out]   lhsValue    The primal value of the lhs. This value is set to the value
       *                           of the right hand side.
       * @param[out]   lhsIndex    The gradient data of the lhs. The index will be set to zero.
       * @param[in]         rhs    The right hand side expression of the assignment.
       */
      CODI_INLINE void store(Real& lhsValue, Index& lhsIndex, const typename TypeTraits<Real>::PassiveReal& rhs) {
        cast().indexHandler.freeIndex(lhsIndex);
        lhsValue = rhs;
      }

      /**
       * @brief Manual store routine.
       *
       * Use this routine to add a statement if the corresponding jacobi entries will be manually pushed onto the tape.
       *
       * The Jacobi entries must be pushed immediately after calling this routine using pushJacobiManual.
       *
       * See also the documentation in TapeInterfaceReverse::storeManual.
       *
       * @param[in]    lhsValue  The primal value of the lhs.
       * @param[out]   lhsIndex  The gradient data of the lhs.
       * @param[in]        size  The number of Jacobi entries.
       */
      CODI_INLINE void storeManual(const Real& lhsValue, Index& lhsIndex, StatementInt size) {
        CODI_UNUSED(lhsValue);

        Tape& tape = cast();

        stmtVector.reserveItems(1);
        tape.jacobiVector.reserveItems(size);
        tape.indexHandler.assignIndex(lhsIndex);
        tape.pushStmtData(size, lhsIndex);
      }

      /**
       * @brief Set the primal value in the primal value vector.
       *
       * Unused in this tape implementation.
       *
       * @param[in]  index  Unused
       * @param[in] primal  Unused
       */
      void setPrimalValue(const Index& index, const Real& primal) {
        CODI_UNUSED(index);
        CODI_UNUSED(primal);
      }


      /**
       * @brief Return the number of used statements.
       * @return The number of used statements.
       */
      size_t getUsedStatementsSize() const {
        return stmtVector.getDataSize();
      }
  };
}
