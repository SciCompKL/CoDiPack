/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2017 Chair for Scientific Computing (SciComp), TU Kaiserslautern
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

/*
 * In order to include this file the user has to define the preprocessor macro CHILD_VECTOR_TYPE,
 * JACOBI_VECTOR_NAME VECTOR_TYPE and STATEMENT_PUSH_FUNCTION_NAME.
 *
 * CHILD_VECTOR_TYPE defines the type of the nested vector for the data vector.
 * JACOBI_VECTOR_NAME defines the member name of the jacobi vector.
 * VECTOR_TYPE defines the type of the data vector.
 * STATEMENT_PUSH_FUNCTION_NAME defines the name of the function that pushes the statements.
 *
 * All these macros are undefined at the end of the file.
 *
 * The module defines the structures stmtVector.
 * The module defines the types StmtChildVector, StmtChildPosition, StmtVector, StmtChunk,
 * StmtPosition.
 *
 * It defines the methods store(Expr), store(const), store(User), printStmtStatistics from the TapeInterface and ReverseTapeInterface.
 *
 * It defines the methods setStatementChunkSize, getUsedStatementSize, resizeStmt as interface functions for the
 * including class.
 */

#ifndef CHILD_VECTOR_TYPE
  #error Please define the type of the child vector
#endif
#ifndef JACOBI_VECTOR_NAME
  #error Please define the name of the Jacobi vector
#endif
#ifndef VECTOR_TYPE
  #error Please define the name of the chunk vector type.
#endif
#ifndef STATEMENT_PUSH_FUNCTION_NAME
  #error Please define the name of the statement push function
#endif


  private:

  // ----------------------------------------------------------------------
  // All definitions of the module
  // ----------------------------------------------------------------------

    /** @brief The child vector for the statement data vector. */
    typedef CHILD_VECTOR_TYPE StmtChildVector;

    /** @brief The position type of the statement child vector */
    typedef typename StmtChildVector::Position StmtChildPosition;

    /** @brief The chunk vector for the statement data. */
    typedef VECTOR_TYPE StmtVector;
    /** @brief The data for each statement. */
    typedef typename StmtVector::ChunkType StmtChunk;

    /** @brief The position type of the statement module. */
    typedef typename StmtVector::Position StmtPosition;

    /** @brief The data for the statements. */
    StmtVector stmtVector;

  private:

  // ----------------------------------------------------------------------
  // Private function for the communication with the including class
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
    CODI_INLINE void store(Real& lhsValue, IndexType& lhsIndex, const Rhs& rhs) {
      void* null = NULL;
      ENABLE_CHECK (OptTapeActivity, active){
        stmtVector.reserveItems(1);
        JACOBI_VECTOR_NAME.reserveItems(ExpressionTraits<Rhs>::maxActiveVariables);
        /* first store the size of the current stack position and evaluate the
         rhs expression. If there was an active variable on the rhs, update
         the index of the lhs */
        size_t startSize = JACOBI_VECTOR_NAME.getChunkPosition();
        rhs.template calcGradient<void*>(null);
        rhs.template pushLazyJacobies<void*>(null);
        size_t activeVariables = JACOBI_VECTOR_NAME.getChunkPosition() - startSize;
        ENABLE_CHECK(OptCheckEmptyStatements, 0 != activeVariables) {

          indexHandler.assignIndex(lhsIndex);
          STATEMENT_PUSH_FUNCTION_NAME((StatementInt)activeVariables, lhsIndex);

#if CODI_AdjointHandle
          Real* jacobies = NULL;
          IndexType* rhsIndices = NULL;

          auto pos = JACOBI_VECTOR_NAME.getPosition();
          JACOBI_VECTOR_NAME.getDataAtPosition(pos.chunk, startSize, jacobies, rhsIndices);

          handleAdjointOperation(rhs.getValue(), lhsIndex, jacobies, rhsIndices, activeVariables);
#endif
        } else {
          indexHandler.freeIndex(lhsIndex);
        }
      } else {
        indexHandler.freeIndex(lhsIndex);
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
    CODI_INLINE void store(Real& lhsValue, IndexType& lhsIndex, const typename TypeTraits<Real>::PassiveReal& rhs) {
      indexHandler.freeIndex(lhsIndex);
      lhsValue = rhs;
    }

    /**
     * @brief Manual store routine.
     *
     * Use this routine to add a statement if the corresponding jacobi entries will be manually pushed onto the tape.
     *
     * The Jacobi entries must be pushed immediately after calling this routine using pushJacobi.
     *
     * @param[in]    lhsValue  The primal value of the lhs.
     * @param[out]   lhsIndex  The gradient data of the lhs.
     * @param[in]        size  The number of Jacobi entries.
     */
    CODI_INLINE void store(const Real& lhsValue, IndexType& lhsIndex, StatementInt size) {
      CODI_UNUSED(lhsValue);

      ENABLE_CHECK (OptTapeActivity, active){
        stmtVector.reserveItems(1);
        JACOBI_VECTOR_NAME.reserveItems(size);
        indexHandler.assignIndex(lhsIndex);
        STATEMENT_PUSH_FUNCTION_NAME(size, lhsIndex);
      }
    }

    /**
     * @brief Prints statistics about the statements.
     *
     * Displays the number of chunks, the total number of statements, the
     * allocated memory and the used memory.
     *
     * @param[in,out]   out  The information is written to the stream.
     * @param[in]     hLine  The horizontal line that separates the sections of the output.
     *
     * @tparam Stream The type of the stream.
     */
    template<typename Stream>
    void printStmtStatistics(Stream& out, const std::string hLine) const {
      size_t nChunksStmts  = stmtVector.getNumChunks();
      size_t totalStmts    = stmtVector.getDataSize();
      size_t sizeStmtEntry = StmtChunk::EntrySize;

      double  memoryUsedStmts = (double)totalStmts*(double)sizeStmtEntry* BYTE_TO_MB;
      double  memoryAllocStmts= (double)nChunksStmts*(double)stmtVector.getChunkSize()
                                *(double)sizeStmtEntry* BYTE_TO_MB;
      out << hLine
          << "Statements\n"
          << hLine
          << "  Total Number:     " << std::setw(10) << totalStmts   << "\n"
          << "  Number of Chunks: " << std::setw(10) << nChunksStmts << "\n"
          << "  Memory used:      " << std::setiosflags(std::ios::fixed)
                                    << std::setprecision(2)
                                    << std::setw(10)
                                    << memoryUsedStmts << " MB" << "\n"
          << "  Memory allocated: " << std::setiosflags(std::ios::fixed)
                                    << std::setprecision(2)
                                    << std::setw(10)
                                    << memoryAllocStmts << " MB" << "\n";

    }

    /**
     * @brief Return the number of used statements.
     * @return The number of used statements.
     */
    size_t getUsedStatementsSize() const {
      return stmtVector.getDataSize();
    }


#undef CHILD_VECTOR_TYPE
#undef JACOBI_VECTOR_NAME
#undef VECTOR_TYPE
#undef STATEMENT_PUSH_FUNCTION_NAME
