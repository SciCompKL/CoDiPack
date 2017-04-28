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
 * The module defines the structures stmtVector, indexVector, constantValueVector, primals, primalsSize, primalsIncr.
 * The module defines the static structures InputHandle, CopyHandle and PreaccHandles,
 * The module defines the types PrimalChildVector, PrimalChildPosition, StatmentVector, StatementChunk, StmtPosition, IndexVector, IndexChunk,
 * IndexPosition, ConstantValueVector, ConstantValueChunk, ConstantValuePosition.
 *
 * It defines the methods store(Expr), store(const), store(User), pushJacobi, printPrimalValueStatistics from the TapeInterface and ReverseTapeInterface.
 *
 * It defines the methods resizePrimals, checkPrimalsSize, evaluateHandle, evaluateConstantValues, getUsedStatementsSize,
 * getUsedDataEntiresSize, getUsedConstantDataSize, setConstantDataSize, swapPrimalValueModule as interface functions for the including class.
 *
 * It defines the static methods inputHandleFunc, copyHandleFunc, preaccHandleFunc as interface functions for the tape.
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

    /** @brief The vector for the statement data. */
    typedef STMT_VECTOR_TYPE StatementVector;
    /** @brief The data for the statments */
    typedef typename StatementVector::ChunkType StatementChunk;
    /** @brief The position type of the statment data. */
    typedef typename StatementVector::Position StmtPosition;
    /** @brief The data for the each statements. */
    StatementVector stmtVector;

    /** @brief The vector for the index data. */
    typedef INDEX_VECTOR_TYPE IndexVector;
    /** @brief The data for the indices*/
    typedef typename IndexVector::ChunkType IndexChunk;
    /** @brief The position type of the index data. */
    typedef typename IndexVector::Position IndexPosition;
    /** @brief The data for the indices of each statements. */
    IndexVector indexVector;

    /** @brief The vector for the constant data. */
    typedef CONSTANT_VECTOR_TYPE ConstantValueVector;
    /** @brief The data for the constants*/
    typedef typename ConstantValueVector::ChunkType ConstantValueChunk;
    /** @brief The position type of the constant data. */
    typedef typename ConstantValueVector::Position ConstantValuePosition;
    /** @brief The data for the constant values of each statements. */
    ConstantValueVector constantValueVector;



    /** @brief The static array of handles for the preaccumulation. */
    const static typename TapeTypes::HandleType preaccHandles[MaxStatementIntSize];

    /** @brief The vector with the primal value for each statement. */
    Real* primals;

    /** @brief The current size of the primal vector. */
    IndexType primalsSize;

    /**
     * @brief The increment of the size of the primal vector.
     *
     * The primals are incremented, when the indices are bigger than the current size.
     */
    IndexType primalsIncr;

  private:

  // ----------------------------------------------------------------------
  // Private function for the communication with the including class
  // ----------------------------------------------------------------------

    /**
     * @brief Swap the data of the primal value module with the data of the other primal tape module.
     *
     * @param[in] other  The object with the other primal value module.
     */
    void swapPrimalValueModule(TAPE_NAME<TapeTypes>& other) {
      std::swap(primals, other.primals);
      std::swap(primalsSize, other.primalsSize);
      std::swap(primalsIncr, other.primalsIncr);
    }

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

    /**
     * @brief Helper function: Cheks if the primal vector is big enough. Increases the vector size otherwise.
     *
     * @param[in] size The new size for the primal vector.
     */
    CODI_INLINE void checkPrimalsSize() {
      if(primalsSize <= indexHandler.getMaximumGlobalIndex()) {
        IndexType newSize = 1 + (indexHandler.getMaximumGlobalIndex() + 1) / primalsIncr;
        newSize = newSize * primalsIncr;
        resizePrimals(newSize);
      }
    }

    /**
     * @brief Action that pushes all passive values to the constant value vector.
     *
     * The function is used on the expression templates as an action.
     *
     * @param[in]  data  Not used in the action.
     * @param[in] value  The passive value that is pushed on the vector.
     */
    CODI_INLINE void pushPassive(int data, const PassiveReal& value) {
      CODI_UNUSED(data);
      constantValueVector.setDataAndMove(value);
    }

    /**
     * @brief Action that counts all active values, that is they have a non zero index.
     *
     * @param[in,out] count  The actual count of the active variables.
     * @param[in]     value  Not used in the action.
     * @param[in]     index  The index of the active real.
     */
    CODI_INLINE void countActiveValues(int* count, const Real& value, const IndexType& index) {
      CODI_UNUSED(value);
      if(0 != index) {
        *count += 1;
      }
    }

    /**
     * @brief Action that adds all indices to the index vector.
     *
     * For passive indices one of the temporary indices is used. e.g 1 to 255.
     * The correspoinding primal value is then pushed to the constant value vector.
     * The constant value is then used in the reverse sweep.
     *
     * @param[in,out] count  The current count of the  inactive variables.
     * @param[in]     value  The primal value of the active real.
     * @param[in]     index  The index of the active real.
     */
    CODI_INLINE void pushIndices(int* passiveVariableCount, const Real& value, const IndexType& index) {
      IndexType pushIndex = index;
      if(0 == pushIndex) {
        *passiveVariableCount += 1;
        pushIndex = *passiveVariableCount;
        constantValueVector.setDataAndMove(value);
      }

      indexVector.setDataAndMove(pushIndex);
    }

  public:
    /**
     * @brief Evaluate one handle in the reverse sweep.
     *
     * The function sets the primal values in the primal value vector for the inactive values.
     * Then it updates the genaral positions and calls the adjoint function of the handle.
     *
     * @param[in]          funcObj  The function object that performs the reverse AD evaluation of an expression.
     * @param[in]          varSize  The number of variables of the expression.
     * @param[in]        constSize  The number constant variables of the expression.
     * @param[in]              adj  The seed from the lhs of the statement.
     * @param[in]   passiveActives  The number of inactive values in the statement.
     * @param[in,out]     indexPos  The position in the index array.
     * @param[in]          indices  The index array.
     * @param[in,out]  constantPos  The position in the constant value array.
     * @param[in]        constants  The constant value array.
     * @param[in,out] primalVector  The global vector with the primal variables.
     * @param[in,out]     adjoints  The adjoint vector for the reverse AD evaluation.
     */
    template<typename FuncObj>
    static CODI_INLINE void evaluateHandle(FuncObj funcObj, size_t varSize, size_t constSize, const GradientValue& adj, const StatementInt& passiveActives, size_t& indexPos, IndexType* &indices, size_t& constantPos, PassiveReal* &constants, Real* primalVector, Real* adjoints) {
      // first restore the primal values of the passive indices
      constantPos -= passiveActives;
      for(StatementInt i = 0; i < passiveActives; ++i) {
        primalVector[i + 1] = constants[constantPos + i];
      }

      // now update the regular pointers
      indexPos -= varSize;
      constantPos -= constSize;
      ENABLE_CHECK(OptZeroAdjoint, adj != 0.0){
        funcObj(adj, &indices[indexPos], &constants[constantPos], primalVector, adjoints);
      }
    }

    private:

    /**
     * @brief Evaluate a part of the index vector.
     *
     * It has to hold start >= end.
     *
     * The function calls the evaluation method for the statement vector.
     *
     * @param[in]    start  The starting point for the index vector.
     * @param[in]      end  The ending point for the index vector.
     * @param[in,out] args  Other arguments for the following functions.
     *
     * @tparam Args  The types of the other arguments.
     */
    template<typename ... Args>
    CODI_INLINE void evaluateIndices(const IndexPosition& start, const IndexPosition& end, Args&&... args) {
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
     * @brief Evaluate a part of the constant value vector.
     *
     * It has to hold start >= end.
     *
     * The function calls the evaluation method for the index vector.
     *
     * @param[in]    start  The starting point for the constant value vector.
     * @param[in]      end  The ending point for the constant value vector.
     * @param[in,out] args  Other arguments for the following functions.
     *
     * @tparam Args  The types of the other arguments.
     */
    template<typename ... Args>
    CODI_INLINE void evaluateConstantValues(const ConstantValuePosition& start, const ConstantValuePosition& end, Args&&... args) {
      PassiveReal* data;
      size_t dataPos = start.data;
      auto curInnerPos = start.inner;
      for(size_t curChunk = start.chunk; curChunk > end.chunk; --curChunk) {
        constantValueVector.getDataAtPosition(curChunk, 0, data);

        auto endInnerPos = constantValueVector.getInnerPosition(curChunk);
        evaluateIndices(curInnerPos, endInnerPos, dataPos, data, std::forward<Args>(args)...);

        codiAssert(dataPos == 0); // after a full chunk is evaluated, the data position needs to be zero

        curInnerPos = endInnerPos;

        dataPos = constantValueVector.getChunkUsedData(curChunk - 1);
      }

      // Iterate over the reminder also covers the case if the start chunk and end chunk are the same
      constantValueVector.getDataAtPosition(end.chunk, 0, data);
      evaluateIndices(curInnerPos, end.inner, dataPos, data, std::forward<Args>(args)...);
    }

    /**
     * @brief Push a copy handle to the tape and all data that is rquired by the handle.
     *
     * The index of the lhs is updated.
     *
     * @param[in]     lhsValue  The value of the active real.
     * @param[in,out] lhsIndex  The index of the active real. The index is renewed in this method.
     * @param[in]     rhsIndex  The index of the rhs value.
     */
    CODI_INLINE void pushCopyHandle(const Real& lhsValue, IndexType& lhsIndex, const IndexType& rhsIndex) {
      indexVector.reserveItems(1);
      indexVector.setDataAndMove(rhsIndex);

      pushStmtData(lhsIndex, lhsValue, HandleFactory::template createHandle<CopyExpr<Real>, TAPE_NAME<TapeTypes> >(), StatementInt(0));
    }

  public:

  // ----------------------------------------------------------------------
  // Public function from the TapeInterface and ReverseTapeInterface
  // ----------------------------------------------------------------------

    /**
     * @brief Store the indices of the statement on the tape.
     *
     * The indices of the rhs expression are stored on the tape.
     * In addition the constant values of the expression are stored in
     * the constant vector and the passive arguments are also stored in the
     * constant vector. For the statement the handle of the rhs expression
     * is written to the tape.
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
          rhs.constantValueAction(*this, NULL, &TAPE_NAME<TapeTypes>::pushPassive);
          codiAssert(ExpressionTraits<Rhs>::maxConstantVariables == constantValueVector.getChunkPosition() - constantSize);

          indexVector.reserveItems(ExpressionTraits<Rhs>::maxActiveVariables);
          size_t indexSize = indexVector.getChunkPosition();
          int passieveVariableCount = 0;
          rhs.valueAction(&passieveVariableCount, &TAPE_NAME<TapeTypes>::pushIndices);
          codiAssert(ExpressionTraits<Rhs>::maxActiveVariables == indexVector.getChunkPosition() - indexSize);
          codiAssert(passieveVariableCount == passiveVariableNumber);

          pushStmtData(lhsIndex, rhs.getValue(), HandleFactory::template createHandle<Rhs, TAPE_NAME<TapeTypes> >(), passiveVariableNumber);

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
     * @param[out]   lhsValue  The primal value of the lhs. This value is set to the value
     *                         of the right hand side.
     * @param[out]   lhsIndex  The gradient data of the lhs. The index will be set to zero.
     * @param[in]         rhs  The right hand side expression of the assignment.
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
     * @param[out]   lhsValue  The primal value of the lhs.
     * @param[out]   lhsIndex  The gradient data of the lhs.
     * @param[in]        size  The number of Jacobi entries.
     */
    CODI_INLINE void store(const Real& lhsValue, IndexType& lhsIndex, StatementInt size) {
      ENABLE_CHECK (OptTapeActivity, active){
        constantValueVector.reserveItems(size);
        indexVector.reserveItems(size);

        pushStmtData(lhsIndex, lhsValue, preaccHandles[size], 0);
      }
    }

    /**
     * @brief Not used in this implementation.
     *
     * @param[in]  data  Not used in this implementation.
     * @param[in] value  Not used in this implementation.
     * @param[in] index  Not used in this implementation.
     */
    template<typename Data>
    CODI_INLINE void pushJacobi(Data& data, const Real& value, const IndexType& index) {
      CODI_UNUSED(data);
      CODI_UNUSED(value);
      CODI_UNUSED(index);

      codiAssert(false || "Should not be called.");
    }

    /**
     * @brief Only used in combination with the manual store.
     *
     * The indices are stored in the index vector. The jacobies
     * are stored in the constant vector.
     *
     * @param[in]     data  Not used in this implementation.
     * @param[in]   jacobi  Stored in the constant vector.
     * @param[in]    value  Not used in this implementation.
     * @param[in]    index  Stored in the index vector.
     */
    template<typename Data>
    CODI_INLINE void pushJacobi(Data& data, const Real& jacobi, const Real& value, const IndexType& index) {
      CODI_UNUSED(data);
      CODI_UNUSED(value);
      CODI_UNUSED(index);

      constantValueVector.setDataAndMove(jacobi);
      indexVector.setDataAndMove(index);
    }

    /**
     * @brief Prints statistics about the stored data for the primal value tape.
     *
     * Displays the number of chunks, the total number of entries, the
     * allocated memory and the used memory for the constant vector, the statement
     * vector and the index vector.
     *
     * @param[in,out]   out  The information is written to the stream.
     * @param[in]     hLine  The horizontal line that separates the sections of the output.
     *
     * @tparam Stream  The type of the stream.
     */
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
          << "Primal Vector\n"
          << hLine
          << "  Total Number:     " << std::setw(10) << totalPrimal << "\n"
          << "  Memory allocated: " << std::setiosflags(std::ios::fixed)
                                    << std::setprecision(2)
                                    << std::setw(10)
                                    << memoryAllocPrimal << " MB" << "\n";
      out << hLine
          << "Statements\n"
          << hLine
          << "  Total Number:     " << std::setw(10) << totalStmt   << "\n"
          << "  Number of Chunks: " << std::setw(10) << nChunksStmt << "\n"
          << "  Memory used:      " << std::setiosflags(std::ios::fixed)
                                    << std::setprecision(2)
                                    << std::setw(10)
                                    << memoryUsedStmt << " MB" << "\n"
          << "  Memory allocated: " << std::setiosflags(std::ios::fixed)
                                    << std::setprecision(2)
                                    << std::setw(10)
                                    << memoryAllocStmt << " MB" << "\n";
      out << hLine
          << "Index entries\n"
          << hLine
          << "  Total Number:     " << std::setw(10) << totalIndex << "\n"
          << "  Number of Chunks: " << std::setw(10) << nChunksIndex << "\n"
          << "  Memory used:      " << std::setiosflags(std::ios::fixed)
                                    << std::setprecision(2)
                                    << std::setw(10)
                                    << memoryUsedIndex << " MB" << "\n"
          << "  Memory allocated: " << std::setiosflags(std::ios::fixed)
                                    << std::setprecision(2)
                                    << std::setw(10)
                                    << memoryAllocIndex << " MB" << "\n";
      out << hLine
          << "Passive data entries\n"
          << hLine
          << "  Total Number:     " << std::setw(10) << totalPassive << "\n"
          << "  Number of Chunks: " << std::setw(10) << nChunksPassive << "\n"
          << "  Memory used:      " << std::setiosflags(std::ios::fixed)
                                    << std::setprecision(2)
                                    << std::setw(10)
                                    << memoryUsedPassive << " MB" << "\n"
          << "  Memory allocated: " << std::setiosflags(std::ios::fixed)
                                    << std::setprecision(2)
                                    << std::setw(10)
                                    << memoryAllocPassive << " MB" << "\n";
    }

    /**
     * @brief Return the number of used stement entries
     * @return The number of used statement entries.
     */
    size_t getUsedStatementsSize() const {
      return stmtVector.getDataSize();
    }

    /**
     * @brief Return the number of used index entries
     * @return The number of used index entries.
     */
    size_t getUsedDataEntriesSize() const {
      return indexVector.getDataSize();
    }

    /**
     * @brief Return the number of used constant data entries
     * @return The number of used constant data entries.
     */
    size_t getUsedConstantDataSize() const {
      return constantValueVector.getDataSize();
    }

    /**
     * @brief Set the number of constant data etnres that should be available.
     *
     * @param[in] constantDataSize  The number of entries that should be available.
     */
    void setConstantDataSize(const size_t& constantDataSize) {
      constantValueVector.resize(constantDataSize);
    }

#undef CHILD_VECTOR_TYPE
#undef VECTOR_TYPE
