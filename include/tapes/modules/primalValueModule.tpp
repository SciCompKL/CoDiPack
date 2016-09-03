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

/*
 * In order to include this file the user has to define the preprocessor macro CHILD_VECTOR_TYPE, CONSTANT_VECTOR_TYPE,
 * INDEX_VECTOR_TYPE, STMT_VECTOR_TYPE and TAPE_NAME.
 *
 * CHILD_VECTOR_TYPE defines the type of the nested vector for the data vector.
 * CONSTANT_VECTOR_TYPE defines the type of the constant value vector.
 * INDEX_VECTOR_TYPE defines the type of the index vector.
 * STMT_TYPE defines the type of the statement vector.
 * TAPE_NAME defines the name of the including tape.
 *
 * All these macros, except TAPE_NAME, are undefined at the end of the file.
 *
 * The module defines the structures TODO.
 * The module defines the types TODO.
 *
 * It defines the methods TODO from the TapeInterface and ReverseTapeInterface.
 *
 * It defines the methods TODO  as interface functions for the
 * including class.
 */

#ifndef CHILD_VECTOR_TYPE
  #error Please define the type of the child vector
#endif
#ifndef CONSTANT_VECTOR_TYPE
  #error Please define the name of the constant value vector type.
#endif
#ifndef INDEX_VECTOR_TYPE
  #error Please define the name of the index vector type.
#endif
#ifndef STMT_VECTOR_TYPE
  #error Please define the name of the statment vector type.
#endif
#ifndef TAPE_NAME
  #error Please define the name of the tape.
#endif


		private:

  // ----------------------------------------------------------------------
  // All definitions of the module
  // ----------------------------------------------------------------------

    /** @brief The child vector for the primal data vector. */
    typedef CHILD_VECTOR_TYPE PrimalChildVector;

    /** @brief The position type of the primal child vector */
    typedef typename PrimalChildVector::Position PrimalChildPosition;

    /** @brief The vector for the primal statement data. */
    typedef STMT_VECTOR_TYPE StatementVector;
    /** @brief The vector for the primal constant data. */
    typedef CONSTANT_VECTOR_TYPE ConstantValueVector;
    /** @brief The vector for the primal index data. */
    typedef INDEX_VECTOR_TYPE IndexVector;

    /** @brief The data for the statments */
    typedef typename StatementVector::ChunkType StatementChunk;
    /** @brief The data for the constants*/
    typedef typename ConstantValueVector::ChunkType ConstantValueChunk;
    /** @brief The data for the indices*/
    typedef typename IndexVector::ChunkType IndexChunk;

    /** @brief The position type of the statment data. */
    typedef typename StatementVector::Position StmtPosition;
    /** @brief The position type of the constant data. */
    typedef typename ConstantValueVector::Position ConstantValuePosition;
    /** @brief The position type of the index data. */
    typedef typename IndexVector::Position IndexPosition;

    /** @brief The data for the each statements. */
    StatementVector stmtVector;
    /** @brief The data for the indices of each statements. */
    IndexVector indexVector;
    /** @brief The data for the constant values of each statements. */
    ConstantValueVector constantValueVector;


    static void inputHandleFunc(const Real& seed, const IndexType* indices, const PassiveReal* constantValues, const Real* primalValues, Real* adjointValues) {}
    const static ExpressionHandle<Real*, Real, IndexType> InputHandle;

    static void copyHandleFunc(const Real& seed, const IndexType* indices, const PassiveReal* constantValues, const Real* primalValues, Real* adjointValues) {
      CODI_UNUSED(constantValues);
      CODI_UNUSED(primalValues);
      adjointValues[indices[0]] += seed;
    }
    const static ExpressionHandle<Real*, Real, IndexType> CopyHandle;

    template<int size>
    static void preaccHandleFunc(const Real& seed, const IndexType* indices, const PassiveReal* constantValues, const Real* primalValues, Real* adjointValues) {
      CODI_UNUSED(primalValues);
      for(int i = 0; i < size; ++i) {
        // jacobies are stored in the constant values
        adjointValues[indices[i]] += constantValues[i] * seed;
      }
    }
    const static ExpressionHandle<Real*, Real, IndexType> PreaccHandles[MaxStatementIntSize];

    Real* primals;
    IndexType primalsSize;
    IndexType primalsIncr;

  private:

  // ----------------------------------------------------------------------
  // Private function for the communication with the including class
  // ----------------------------------------------------------------------

    /**
     * @brief Helper function: Sets the primal vector to a new size.
     *
     * @param[in] size The new size for the primal vector.
     */
    void resizePrimals(const IndexType& size) {
      IndexType oldSize = primalsSize;
      primalsSize = size;

      primals = (Real*)realloc(primals, sizeof(Real) * (size_t)primalsSize);

      for(IndexType i = oldSize; i < primalsSize; ++i) {
        primals[i] = Real();
      }
    }

    CODI_INLINE void checkPrimalsSize() {
      if(primalsSize <= indexHandler.getMaximumGlobalIndex()) {
        IndexType newSize = 1 + (indexHandler.getMaximumGlobalIndex() + 1) / primalsIncr;
        newSize = newSize * primalsIncr;
        resizePrimals(newSize);
      }
    }

    CODI_INLINE void pushPassive(int data, const PassiveReal& value) {
      CODI_UNUSED(data);
      constantValueVector.setDataAndMove(value);
    }

    CODI_INLINE void countActiveValues(int* count, const Real& value, const IndexType& index) {
      CODI_UNUSED(value);
      if(0 != index) {
        *count += 1;
      }
    }

    CODI_INLINE void pushIndices(int* passiveVariableCount, const Real& value, const IndexType& index) {
      IndexType pushIndex = index;
      if(0 == pushIndex) {
        *passiveVariableCount += 1;
        pushIndex = *passiveVariableCount;
        constantValueVector.setDataAndMove(value);
      }

      indexVector.setDataAndMove(pushIndex);
    }

    CODI_INLINE void evaluateHandle(const GradientValue& adj, Handle& exprHandle, StatementInt& passiveActives, size_t& indexPos, IndexType* &indices, size_t& constantPos, PassiveReal* &constants) {
      // first restore the primal values of the passive indices
      constantPos -= passiveActives;
      for(StatementInt i = 0; i < passiveActives; ++i) {
        primals[i + 1] = constants[constantPos + i];
      }

      // now update the regular pointers
      indexPos -= exprHandle->maxActiveVariables;
      constantPos -= exprHandle->maxConstantVariables;
      ENABLE_CHECK(OptZeroAdjoint, adj != 0.0){

        exprHandle->adjointFunc(adj, &indices[indexPos], &constants[constantPos], primals, adjoints);
      }
    }

    template<typename ... Args>
    CODI_INLINE void evalIndices(const IndexPosition& start, const IndexPosition& end, Args&&... args) {
      IndexType* data;
      size_t dataPos = start.data;
      auto curInnerPos = start.inner;
      for(size_t curChunk = start.chunk; curChunk > end.chunk; --curChunk) {
        indexVector.getDataAtPosition(curChunk, 0, data);

        auto endInnerPos = indexVector.getInnerPosition(curChunk);
        evalStmt(curInnerPos, endInnerPos, dataPos, data, std::forward<Args>(args)...);

        codiAssert(dataPos == 0); // after a full chunk is evaluated, the data position needs to be zero

        curInnerPos = endInnerPos;

        dataPos = indexVector.getChunkUsedData(curChunk - 1);
      }

      // Iterate over the reminder also covers the case if the start chunk and end chunk are the same
      indexVector.getDataAtPosition(end.chunk, 0, data);
      evalStmt(curInnerPos, end.inner, dataPos, data, std::forward<Args>(args)...);
    }

    /**
     * @brief Evaluate a part of the statement vector.
     *
     * It has to hold start >= end.
     *
     * The function calls the evaluation method for the jacobi vector.
     *
     * @param[in] start The starting point for the statement vector.
     * @param[in]   end The ending point for the statement vector.
     */
    template<typename ... Args>
    CODI_INLINE void evalExtFuncCallback(const ConstantValuePosition& start, const ConstantValuePosition& end, Args&&... args) {
      PassiveReal* data;
      size_t dataPos = start.data;
      auto curInnerPos = start.inner;
      for(size_t curChunk = start.chunk; curChunk > end.chunk; --curChunk) {
        constantValueVector.getDataAtPosition(curChunk, 0, data);

        auto endInnerPos = constantValueVector.getInnerPosition(curChunk);
        evalIndices(curInnerPos, endInnerPos, dataPos, data, std::forward<Args>(args)...);

        codiAssert(dataPos == 0); // after a full chunk is evaluated, the data position needs to be zero

        curInnerPos = endInnerPos;

        dataPos = constantValueVector.getChunkUsedData(curChunk - 1);
      }

      // Iterate over the reminder also covers the case if the start chunk and end chunk are the same
      constantValueVector.getDataAtPosition(end.chunk, 0, data);
      evalIndices(curInnerPos, end.inner, dataPos, data, std::forward<Args>(args)...);
    }

    CODI_INLINE void pushCopyHandle(const Real& lhsValue, IndexType& lhsIndex, const IndexType& rhsIndex) {
      indexVector.reserveItems(1);
      indexVector.setDataAndMove(rhsIndex);

      pushStmtData(lhsIndex, lhsValue, &CopyHandle, StatementInt(0));
    }

  public:

  // ----------------------------------------------------------------------
  // Public function from the TapeInterface and ReverseTapeInterface
  // ----------------------------------------------------------------------

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

      ENABLE_CHECK(OptTapeActivity, active){

        int activeCount = 0;
        rhs.valueAction(&activeCount, &TAPE_NAME<TapeTypes>::countActiveValues);

        if(0 != activeCount) {
          int passiveVariableNumber = ExpressionTraits<Rhs>::maxActiveVariables - activeCount;

          constantValueVector.reserveItems(ExpressionTraits<Rhs>::maxConstantVariables + passiveVariableNumber); // the additional passives are create in pushIndices
          size_t constantSize = constantValueVector.getChunkPosition();
          rhs.passiveAction(*this, NULL, &TAPE_NAME<TapeTypes>::pushPassive);
          codiAssert(ExpressionTraits<Rhs>::maxConstantVariables == constantValueVector.getChunkPosition() - constantSize);

          indexVector.reserveItems(ExpressionTraits<Rhs>::maxActiveVariables);
          size_t indexSize = indexVector.getChunkPosition();
          int passieveVariableCount = 0;
          rhs.valueAction(&passieveVariableCount, &TAPE_NAME<TapeTypes>::pushIndices);
          codiAssert(ExpressionTraits<Rhs>::maxActiveVariables == indexVector.getChunkPosition() - indexSize);
          codiAssert(passieveVariableCount == passiveVariableNumber);

          pushStmtData(lhsIndex, rhs.getValue(), ExpressionHandleStore<Real*, Real, IndexType, Rhs>::getHandle(), passiveVariableNumber);

          CODI_UNUSED(constantSize);  /* needed to avoid unused variable when the assersts are not enabled. */
          CODI_UNUSED(indexSize);  /* needed to avoid unused variable when the assersts are not enabled. */

#if CODI_AdjointHandle
            IndexType* rhsIndices = NULL;
            PassiveReal* constants = NULL;

            auto posIndex = indexVector.getPosition();
            indexVector.getDataAtPosition(posIndex.chunk, indexSize, rhsIndices);

            auto posPassive = constantValueVector.getPosition();
            constantValueVector.getDataAtPosition(posPassive.chunk, constantSize, constants);

            resizeAdjoints(indexHandler.getMaximumGlobalIndex() + 1);
            handleAdjointOperation(rhs.getValue(), lhsIndex, ExpressionHandleStore<Real*, Real, IndexType, Rhs>::getHandle(), passiveVariableNumber, constants, rhsIndices, primals, adjoints);
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
     * @param[out]   lhsIndex    The gradient data of the lhs.
     * @param[in]        size    The number of Jacobi entries.
     */
    CODI_INLINE void store(const Real& lhsValue, IndexType& lhsIndex, StatementInt size) {
      ENABLE_CHECK (OptTapeActivity, active){
        constantValueVector.reserveItems(size);
        indexVector.reserveItems(size);

        pushStmtData(lhsIndex, lhsValue, &PreaccHandles[size], 0);
      }
    }

    /**
     * @brief Stores the jacobi with the value 1.0 on the tape if the index is active.
     *
     * @param[in] gradient Not used in this implementation.
     * @param[in]    value Not used in this implementation.
     * @param[in]    index Used to check if the variable is active.
     */
    template<typename Data>
    CODI_INLINE void pushJacobi(Data& data, const Real& value, const IndexType& index) {
      CODI_UNUSED(data);
      CODI_UNUSED(value);
      CODI_UNUSED(index);

      codiAssert(false || "Should not be called.");
    }

    /**
     * @brief Stores the jacobi on the tape if the index is active.
     *
     * @param[in] gradient Not used in this implementation.
     * @param[in]   jacobi Stored on the tape if the variable is active.
     * @param[in]    value Not used in this implementation.
     * @param[in]    index Used to check if the variable is active.
     */
    template<typename Data>
    CODI_INLINE void pushJacobi(Data& data, const Real& jacobi, const Real& value, const IndexType& index) {
      CODI_UNUSED(data);
      CODI_UNUSED(value);
      CODI_UNUSED(index);

      constantValueVector.setDataAndMove(jacobi);
      indexVector.setDataAndMove(index);
    }

    template<typename Stream>
    void printPrimalValueStatistics(Stream& out, const std::string& hLine) const {

      size_t nChunksIndex  = indexVector.getNumChunks();
      size_t totalIndex    = indexVector.getDataSize();
      size_t sizeIndexEntry = IndexChunk::EntrySize;
      double memoryUsedIndex = (double)totalIndex*(double)sizeIndexEntry* BYTE_TO_MB;
      double memoryAllocIndex= (double)nChunksIndex*(double)indexVector.getChunkSize()*(double)sizeIndexEntry* BYTE_TO_MB;

      size_t nChunksStmt  = stmtVector.getNumChunks();
      size_t totalStmt    = stmtVector.getDataSize();
      size_t sizeStmtEntry = StatementChunk::EntrySize;
      double memoryUsedStmt = (double)totalStmt*(double)sizeStmtEntry* BYTE_TO_MB;
      double memoryAllocStmt= (double)nChunksStmt*(double)stmtVector.getChunkSize()*(double)sizeStmtEntry* BYTE_TO_MB;

      size_t nChunksPassive  = constantValueVector.getNumChunks();
      size_t totalPassive    = constantValueVector.getDataSize();
      size_t sizePassiveEntry = ConstantValueChunk::EntrySize;
      double memoryUsedPassive = (double)totalPassive*(double)sizePassiveEntry* BYTE_TO_MB;
      double memoryAllocPassive= (double)nChunksPassive*(double)constantValueVector.getChunkSize()*(double)sizePassiveEntry* BYTE_TO_MB;

      size_t totalPrimal   = primalsSize;
      size_t sizePrimalEntry = sizeof(Real);
      double memoryAllocPrimal = (double)totalPrimal*(double)sizePrimalEntry* BYTE_TO_MB;

      out << hLine
          << "Primal Vector \n"
          << hLine
          << "  Total Number:     " << std::setw(10) << totalPrimal << "\n"
          << "  Memory allocated: " << std::setiosflags(std::ios::fixed)
                                    << std::setprecision(2)
                                    << std::setw(10)
                                    << memoryAllocPrimal << " MB" << "\n";
      out << hLine
          << "Statements \n"
          << hLine
          << "  Number of Chunks: " << std::setw(10) << nChunksStmt << "\n"
          << "  Total Number:     " << std::setw(10) << totalStmt   << "\n"
          << "  Memory allocated: " << std::setiosflags(std::ios::fixed)
                                    << std::setprecision(2)
                                    << std::setw(10)
                                    << memoryAllocStmt << " MB" << "\n"
          << "  Memory used:      " << std::setiosflags(std::ios::fixed)
                                    << std::setprecision(2)
                                    << std::setw(10)
                                    << memoryUsedStmt << " MB" << "\n";
      out << hLine
          << "Index entries \n"
          << hLine
          << "  Number of Chunks: " << std::setw(10) << nChunksIndex << "\n"
          << "  Total Number:     " << std::setw(10) << totalIndex << "\n"
          << "  Memory allocated: " << std::setiosflags(std::ios::fixed)
                                    << std::setprecision(2)
                                    << std::setw(10)
                                    << memoryAllocIndex << " MB" << "\n"
          << "  Memory used:      " << std::setiosflags(std::ios::fixed)
                                    << std::setprecision(2)
                                    << std::setw(10)
                                    << memoryUsedIndex << " MB" << "\n";
      out << hLine
          << "Passive data entries \n"
          << hLine
          << "  Number of Chunks: " << std::setw(10) << nChunksPassive << "\n"
          << "  Total Number:     " << std::setw(10) << totalPassive << "\n"
          << "  Memory allocated: " << std::setiosflags(std::ios::fixed)
                                    << std::setprecision(2)
                                    << std::setw(10)
                                    << memoryAllocPassive << " MB" << "\n"
          << "  Memory used:      " << std::setiosflags(std::ios::fixed)
                                    << std::setprecision(2)
                                    << std::setw(10)
                                    << memoryUsedPassive << " MB" << "\n";
    }

    size_t getUsedStatementsSize() const {
      return stmtVector.getDataSize();
    }

    size_t getUsedDataEntriesSize() const {
      return indexVector.getDataSize();
    }

    size_t getUsedConstantDataSize() const {
      return constantValueVector.getDataSize();
    }

    void setConstantDataSize(const size_t& constantDataSize) {
      constantValueVector.resize(constantDataSize);
    }

#undef CHILD_VECTOR_TYPE
#undef VECTOR_TYPE
