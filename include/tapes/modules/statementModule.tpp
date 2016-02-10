/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015 Chair for Scientific Computing (SciComp), TU Kaiserslautern
 * Homepage: http://www.scicomp.uni-kl.de
 * Contact:  Prof. Nicolas R. Gauger (codi@scicomp.uni-kl.de)
 *
 * Lead developers: Max Sagebaum, Tim Albring (SciComp, TU Kaiserslautern)
 *
 * This file is part of CoDiPack (http://www.scicomp.uni-kl.de/software/codi).
 *
 * CoDiPack is free software: you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation, either version 2 of the
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


  public:

    typedef CHILD_VECTOR_TYPE StmtChildVector;
    typedef typename StmtChildVector::Position StmtChildPosition;

    /** @brief The chunk vector for the statement data. */
    typedef VECTOR_TYPE StmtVector;
    /** @brief The data for each statement. */
    typedef typename StmtVector::ChunkType StmtChunk;

    typedef typename StmtVector::Position StmtPosition;

    /** @brief The data for the statements. */
    StmtVector stmtVector;

    /**
     * @brief Set the size of the statement data chunks.
     *
     * @param[in] statementChunkSize The new size for the statement data chunks.
     */
    void setStatementChunkSize(const size_t& statementChunkSize) {
      stmtVector.setChunkSize(statementChunkSize);
    }

    /**
     * @brief Return the number of used statements.
     * @return The number of used statements.
     */
    size_t getUsedStatementsSize() const {
      return stmtVector.getDataSize();
    }

    void resizeStmt(const size_t& statementSize) {
      stmtVector.resize(statementSize);
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
    inline void store(Real& lhsValue, IndexType& lhsIndex, const Rhs& rhs) {
      void* null = NULL;
      ENABLE_CHECK (OptTapeActivity, active){
        stmtVector.reserveItems(1);
        JACOBI_VECTOR_NAME.reserveItems(ExpressionTraits<Rhs>::maxActiveVariables);
        /* first store the size of the current stack position and evaluate the
         rhs expression. If there was an active variable on the rhs, update
         the index of the lhs */
        size_t startSize = JACOBI_VECTOR_NAME.getChunkPosition();
        rhs.template calcGradient<void*>(null);
        size_t activeVariables = JACOBI_VECTOR_NAME.getChunkPosition() - startSize;
        ENABLE_CHECK(OptCheckEmptyStatements, 0 != activeVariables) {

          indexHandler.checkIndex(lhsIndex);
          STATEMENT_PUSH_FUNCTION_NAME((StatementInt)activeVariables, lhsIndex);
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
    inline void store(Real& lhsValue, IndexType& lhsIndex, const typename TypeTraits<Real>::PassiveReal& rhs) {
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
     * @param[out]   lhsIndex    The gradient data of the lhs.
     * @param[in]        size    The number of Jacobi entries.
     */
    inline void store(IndexType& lhsIndex, StatementInt size) {
      ENABLE_CHECK (OptTapeActivity, active){
        stmtVector.reserveItems(1);
        JACOBI_VECTOR_NAME.reserveItems(size);
        indexHandler.checkIndex(lhsIndex);
        stmtVector.setDataAndMove(std::make_tuple(size));
        STATEMENT_PUSH_FUNCTION_NAME(size, lhsIndex);
      }
    }

    /**
     * @brief Prints statistics about the tape on the screen
     *
     * Prints information such as stored statements/adjoints and memory usage on screen.
     */
    void printStmtStatistics() const {
      size_t nChunksStmts  = stmtVector.getNumChunks();
      size_t totalStmts    = (nChunksStmts-1)*stmtVector.getChunkSize()
                             +stmtVector.getChunkUsedData(nChunksStmts-1);
      double  memoryUsedStmts = (double)totalStmts*(double)sizeof(StatementInt)* BYTE_TO_MB;
      double  memoryAllocStmts= (double)nChunksStmts*(double)stmtVector.getChunkSize()
                                *(double)sizeof(StatementInt)* BYTE_TO_MB;
      std::cout << std::endl
                << "Statements " << std::endl
                << "-------------------------------------" << std::endl
                << "  Number of Chunks: " << std::setw(10) << nChunksStmts << std::endl
                << "  Total Number:     " << std::setw(10) << totalStmts   << std::endl
                << "  Memory allocated: " << std::setiosflags(std::ios::fixed)
                                          << std::setprecision(2)
                                          << std::setw(10)
                                          << memoryAllocStmts << " MB" << std::endl
                << "  Memory used:      " << std::setiosflags(std::ios::fixed)
                                          << std::setprecision(2)
                                          << std::setw(10)
                                          << memoryUsedStmts << " MB" << std::endl;

    }

#undef CHILD_VECTOR_TYPE
#undef JACOBI_VECTOR_NAME
#undef VECTOR_TYPE
#undef STATEMENT_PUSH_FUNCTION_NAME
