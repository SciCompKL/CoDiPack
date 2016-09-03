/**
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


#pragma once

#include <cstddef>
#include <tuple>

#include "../activeReal.hpp"
#include "../expressionHandle.hpp"
#include "chunkVector.hpp"
#include "indices/reuseIndexHandler.hpp"
#include "reverseTapeInterface.hpp"
#include "singleChunkVector.hpp"

namespace codi {

  /**
   * @brief Vector defintion for the ChunkPrimalValueIndexTape.
   *
   * The structure defines all vectors as chunk vectors.
   *
   * See PrimalValueIndexTape for details.
   *
   * @tparam Real  The type for the primal values.
   * @tparam IndexHandler  The index handler for the managing of the indices. It has to be a index handler that assumes index reuse.
   * @tparam GradientValue  The type for the adjoint values. (Default: Same as the primal value.)
   */
  template <typename Real, typename IndexHandler, typename GradientValue = Real>
  struct ChunkIndexPrimalValueTapeTypes {
    /** @brief The type for the primal values. */
    typedef Real RealType;
    /** @brief The handler for the indices. */
    typedef IndexHandler IndexHandlerType;
    /** @brief The type for the adjoint values. */
    typedef GradientValue GradientValueType;

    typedef typename IndexHandler::IndexType IndexType;

    typedef const ExpressionHandle<Real*, Real, IndexType>* HandleType;

    /** @brief The data for each statement. */
    typedef Chunk4<IndexType, Real, HandleType, StatementInt> StatementChunk;
    /** @brief The chunk vector for the statement data. */
    typedef ChunkVector<StatementChunk, EmptyChunkVector> StatementVector;

    /** @brief The data for the indices of each statement */
    typedef Chunk1< typename IndexHandler::IndexType> IndexChunk;
    /** @brief The chunk vector for the index data. */
    typedef ChunkVector<IndexChunk, StatementVector> IndexVector;

    /** @brief The data for the constant values of each statement */
    typedef Chunk1< typename TypeTraits<Real>::PassiveReal> ConstantValueChunk;
    /** @brief The chunk vector for the constant data. */
    typedef ChunkVector<ConstantValueChunk, IndexVector> ConstantValueVector;

    /** @brief The data for the external functions. */
    typedef Chunk2<ExternalFunction,typename ConstantValueVector::Position> ExternalFunctionChunk;
    /** @brief The chunk vector for the external  function data. */
    typedef ChunkVector<ExternalFunctionChunk, ConstantValueVector> ExternalFunctionVector;

    /** @brief The position for all the different data vectors. */
    typedef typename ExternalFunctionVector::Position Position;

    constexpr static const char* tapeName = "ChunkPrimalValueIndexTape";

  };

  /**
   * @brief Vector defintion for the SimplePrimalValueIndexTape.
   *
   * The structure defines all vectors as single chunk vectors.
   *
   * See PrimalValueIndexTape for details.
   *
   * @tparam Real  The type for the primal values.
   * @tparam IndexHandler  The index handler for the managing of the indices. It has to be a index handler that assumes index reuse.
   * @tparam GradientValue  The type for the adjoint values. (Default: Same as the primal value.)
   */
  template <typename Real, typename IndexHandler, typename GradientValue = Real>
  struct SimpleIndexPrimalValueTapeTypes {
    /** @brief The type for the primal values. */
    typedef Real RealType;
    /** @brief The handler for the indices. */
    typedef IndexHandler IndexHandlerType;
    /** @brief The type for the adjoint values. */
    typedef GradientValue GradientValueType;

    typedef typename IndexHandler::IndexType IndexType;

    typedef const ExpressionHandle<Real*, Real, IndexType>* HandleType;

    /** @brief The data for each statement. */
    typedef Chunk4<IndexType, Real, HandleType, StatementInt> StatementChunk;
    /** @brief The chunk vector for the statement data. */
    typedef SingleChunkVector<StatementChunk, EmptyChunkVector> StatementVector;

    /** @brief The data for the indices of each statement */
    typedef Chunk1< typename IndexHandler::IndexType> IndexChunk;
    /** @brief The chunk vector for the index data. */
    typedef SingleChunkVector<IndexChunk, StatementVector> IndexVector;

    /** @brief The data for the constant values of each statement */
    typedef Chunk1< typename TypeTraits<Real>::PassiveReal> ConstantValueChunk;
    /** @brief The chunk vector for the constant data. */
    typedef SingleChunkVector<ConstantValueChunk, IndexVector> ConstantValueVector;

    /** @brief The data for the external functions. */
    typedef Chunk2<ExternalFunction,typename ConstantValueVector::Position> ExternalFunctionChunk;
    /** @brief The chunk vector for the external  function data. */
    typedef SingleChunkVector<ExternalFunctionChunk, ConstantValueVector> ExternalFunctionVector;

    /** @brief The position for all the different data vectors. */
    typedef typename ExternalFunctionVector::Position Position;

    constexpr static const char* tapeName = "SimplePrimalValueIndexTape";

  };

  template <typename TapeTypes>
  class PrimalValueIndexTape : public ReverseTapeInterface<typename TapeTypes::RealType, typename TapeTypes::IndexHandlerType::IndexType, typename TapeTypes::GradientValueType, PrimalValueIndexTape<TapeTypes>, typename TapeTypes::Position > {
  public:

    /** @brief The type for the primal values. */
    typedef typename TapeTypes::RealType Real;
    /** @brief The type for the adjoint values. */
    typedef typename TapeTypes::GradientValueType GradientValue;
    /** @brief The type for the index handler. */
    typedef typename TapeTypes::IndexHandlerType IndexHandler;

    /** @brief The type for the indices that are used for the identification of the adjoint variables. */
    typedef typename TapeTypes::IndexType IndexType;
    /** @brief The gradient data is just the index type. */
    typedef IndexType GradientData;

    typedef typename TypeTraits<Real>::PassiveReal PassiveReal;
    typedef typename TapeTypes::HandleType Handle;

    typedef typename TapeTypes::ConstantValueVector ConstantValueVector;
    typedef typename TapeTypes::IndexVector IndexVector;
    typedef typename TapeTypes::StatementVector StatementVector;

    typedef typename TapeTypes::ConstantValueVector::Position PassivePosition;
    typedef typename TapeTypes::IndexVector::Position IndexPosition;
    typedef typename TapeTypes::StatementVector::Position StmtPosition;

    /** @brief The counter for the current expression. */
    EmptyChunkVector emptyVector;

    #define TAPE_NAME PrimalValueIndexTape

    #define POSITION_TYPE typename TapeTypes::Position
    #define INDEX_HANDLER_TYPE IndexHandler
    #define RESET_FUNCTION_NAME resetInt
    #define EVALUATE_FUNCTION_NAME evaluateInt
    #include "modules/tapeBaseModule.tpp"

    typename TapeTypes::StatementVector stmtVector;
    typename TapeTypes::IndexVector indexVector;
    typename TapeTypes::ConstantValueVector constantValueVector;

    #define CHILD_VECTOR_TYPE ConstantValueVector
    #define CHILD_VECTOR_NAME constantValueVector
    #define VECTOR_TYPE typename TapeTypes::ExternalFunctionVector
    #include "modules/externalFunctionsModule.tpp"

    #undef TAPE_NAME

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
    Real* primalValueCopy;
    IndexType primalsSize;
    IndexType primalsIncr;

  public:
    /**
     * @brief Creates a tape with the size of zero for the data, statements and external functions.
     */
    PrimalValueIndexTape() :
      emptyVector(),
      /* defined in tapeBaseModule */indexHandler(),
      /* defined in tapeBaseModule */adjoints(NULL),
      /* defined in tapeBaseModule */adjointsSize(0),
      /* defined in tapeBaseModule */active(false),
      stmtVector(DefaultChunkSize, emptyVector),
      indexVector(DefaultChunkSize, stmtVector),
      constantValueVector(DefaultChunkSize, indexVector),
      /* defined in externalFunctionsModule */extFuncVector(1000, constantValueVector),
      primals(NULL),
      primalValueCopy(NULL),
      primalsSize(0),
      primalsIncr(DefaultSmallChunkSize){}

    size_t getUsedStatementsSize() const {
      return stmtVector.getDataSize();
    }

    size_t getUsedDataEntriesSize() const {
      return indexVector.getDataSize();
    }

    size_t getUsedPassiveDataSize() const {
      return constantValueVector.getDataSize();
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

    inline void checkPrimalsSize() {
      if(primalsSize <= indexHandler.getMaximumGlobalIndex()) {
        IndexType newSize = 1 + (indexHandler.getMaximumGlobalIndex() + 1) / primalsIncr;
        newSize = newSize * primalsIncr;
        resizePrimals(newSize);
      }
    }

    /**
     * @brief Sets all adjoint/gradients to zero.
     *
     * It has to hold start >= end.
     *
     * @param[in] start  The starting position for the reset of the vector.
     * @param[in]   end  The ending position for the reset of the vector.
     */
    inline void clearAdjoints(const Position& start, const Position& end){
      CODI_UNUSED(start);
      CODI_UNUSED(end);
    }

    void setConstantDataSize(const size_t& constantDataSize) {
      constantValueVector.resize(constantDataSize);
    }

    /**
     * @brief Set the size of the jacobi and statement data and the adjoint vector.
     * @param[in] dataSize  The new size of the jacobi vector.
     * @param[in] stmtSize  The new size of the statement vector.
     */
    void resize(const size_t& dataSize, const size_t& stmtSize) {
      indexVector.resize(dataSize);
      stmtVector.resize(stmtSize);

      resizePrimals(stmtSize + 1);
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

      ENABLE_CHECK(OptTapeActivity, active){

        int activeCount = 0;
        rhs.valueAction(&activeCount, &PrimalValueIndexTape<TapeTypes>::countActiveValues);

        if(0 != activeCount) {
          int passiveVariableNumber = ExpressionTraits<Rhs>::maxActiveVariables - activeCount;

          constantValueVector.reserveItems(ExpressionTraits<Rhs>::maxConstantVariables + passiveVariableNumber); // the additional passives are create in pushIndices
          size_t constantSize = constantValueVector.getChunkPosition();
          rhs.passiveAction(*this, NULL, &PrimalValueIndexTape<TapeTypes>::pushPassive);
          codiAssert(ExpressionTraits<Rhs>::maxConstantVariables == constantValueVector.getChunkPosition() - constantSize);

          indexVector.reserveItems(ExpressionTraits<Rhs>::maxActiveVariables);
          size_t indexSize = indexVector.getChunkPosition();
          int passieveVariableCount = 0;
          rhs.valueAction(&passieveVariableCount, &PrimalValueIndexTape<TapeTypes>::pushIndices);
          codiAssert(ExpressionTraits<Rhs>::maxActiveVariables == indexVector.getChunkPosition() - indexSize);
          codiAssert(passieveVariableCount == passiveVariableNumber);

          indexHandler.assignIndex(lhsIndex);
          stmtVector.reserveItems(1);
          checkPrimalsSize();
          stmtVector.setDataAndMove(lhsIndex, primals[lhsIndex], ExpressionHandleStore<Real*, Real, IndexType, Rhs>::getHandle(), passiveVariableNumber);

          primals[lhsIndex] = rhs.getValue();

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
     * @brief Optimization for the copy operation just copies the index of the rhs.
     *
     * No data is stored in this method.
     *
     * The primal value of the lhs is set to the primal value of the rhs.
     *
     * @param[out]   lhsValue    The primal value of the lhs. This value is set to the value
     *                           of the right hand side.
     * @param[out]   lhsIndex    The gradient data of the lhs. The index will be set to the index of the rhs.
     * @param[in]         rhs    The right hand side expression of the assignment.
     */
    inline void store(Real& lhsValue, IndexType& lhsIndex, const ActiveReal<PrimalValueIndexTape<TapeTypes> >& rhs) {
      ENABLE_CHECK(OptTapeActivity, active){
        ENABLE_CHECK(OptCheckZeroIndex, 0 != rhs.getGradientData()) {
          indexHandler.copyIndex(lhsIndex, rhs.getGradientData());

          if(IndexHandler::AssignNeedsStatement) {
            pushCopyHandle(rhs.getValue(), lhsIndex, rhs.getGradientData());
          }
        } else {
          indexHandler.freeIndex(lhsIndex);
        }
      } else {
        indexHandler.freeIndex(lhsIndex);
      }
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
    CODI_INLINE void store(const Real& lhsValue, IndexType& lhsIndex, StatementInt size) {
      ENABLE_CHECK (OptTapeActivity, active){
        constantValueVector.reserveItems(size);
        indexVector.reserveItems(size);
        stmtVector.reserveItems(1);
        indexHandler.assignIndex(lhsIndex);
        stmtVector.setDataAndMove(lhsIndex, primals[lhsIndex], &PreaccHandles[size], 0);

        primals[lhsIndex] = lhsValue;
      }
    }

    inline void pushPassive(int data, const PassiveReal& value) {
      CODI_UNUSED(data);
      constantValueVector.setDataAndMove(value);
    }

    inline void countActiveValues(int* count, const Real& value, const IndexType& index) {
      CODI_UNUSED(value);
      if(0 != index) {
        *count += 1;
      }
    }

    inline void pushIndices(int* passiveVariableCount, const Real& value, const IndexType& index) {
      IndexType pushIndex = index;
      if(0 == pushIndex) {
        *passiveVariableCount += 1;
        pushIndex = *passiveVariableCount;
        constantValueVector.setDataAndMove(value);
      }

      indexVector.setDataAndMove(pushIndex);
    }

    /**
     * @brief Stores the jacobi with the value 1.0 on the tape if the index is active.
     *
     * @param[in] gradient Not used in this implementation.
     * @param[in]    value Not used in this implementation.
     * @param[in]    index Used to check if the variable is active.
     */
    template<typename Data>
    inline void pushJacobi(Data& data, const Real& value, const IndexType& index) {
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
    inline void pushJacobi(Data& data, const Real& jacobi, const Real& value, const IndexType& index) {
      CODI_UNUSED(data);
      CODI_UNUSED(value);
      CODI_UNUSED(index);

      constantValueVector.setDataAndMove(jacobi);
      indexVector.setDataAndMove(index);
    }

    /**
     * @brief Get the current position of the tape.
     *
     * The position can be used to reset the tape to that position or to
     * evaluate only parts of the tape.
     * @return The current position of the tape.
     */
    inline Position getPosition() const {
      return getExtFuncPosition();
    }

    /**
     * @brief Get the current position of the tape.
     *
     * The position can be used to reset the tape to that position or to
     * evaluate only parts of the tape.
     * @return The current position of the tape.
     */
    inline Position getZeroPosition() const {
      return getExtFuncZeroPosition();
    }

    /**
     * @brief Reset the tape to the given position.
     *
     * @param[in] pos Reset the state of the tape to the given position.
     */
    inline void resetInt(const Position& pos) {
      resetExtFunc(pos);
    }

  private:
    /**
     * @brief Evaluate the stack from the start to to the end position.
     *
     * It has to hold start >= end.
     *
     * @param[in] start The starting position for the adjoint evaluation.
     * @param[in]   end The ending position for the adjoint evaluation.
     */
    inline void evaluateStack(size_t& stmtPos, const size_t& endStmtPos, IndexType* lhsIndices, Real* storedPrimals, Handle* &statements, StatementInt* &passiveActiveReal, size_t& indexPos, IndexType* &indices, size_t& constantPos, PassiveReal* &constants, Real* primalVector) {
      while(stmtPos > endStmtPos) {
        --stmtPos;
        const IndexType& lhsIndex = lhsIndices[stmtPos];
        Handle exprHandle = statements[stmtPos];

        primalVector[lhsIndex] = storedPrimals[stmtPos];
        const GradientValue adj = adjoints[lhsIndex];
        adjoints[lhsIndex] = GradientValue();

        // first restore the primal values of the passive indices
        StatementInt passiveActives = passiveActiveReal[stmtPos];
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
    inline void evalStmt(const StmtPosition& start, const StmtPosition& end, Args&&... args) {
      IndexType* data1;
      Real* data2;
      Handle* data3;
      StatementInt* data4;

      size_t dataPos = start.data;
      for(size_t curChunk = start.chunk; curChunk > end.chunk; --curChunk) {
        stmtVector.getDataAtPosition(curChunk, 0, data1, data2, data3, data4);

        evaluateStack(dataPos, 0, data1, data2, data3, data4, std::forward<Args>(args)..., primalValueCopy);

        dataPos = stmtVector.getChunkUsedData(curChunk - 1);
      }

      // Iterate over the reminder also covers the case if the start chunk and end chunk are the same
      stmtVector.getDataAtPosition(end.chunk, 0, data1, data2, data3, data4);
      evaluateStack(dataPos, end.data, data1, data2, data3, data4, std::forward<Args>(args)..., primalValueCopy);
    }

    template<typename ... Args>
    inline void evalIndices(const IndexPosition& start, const IndexPosition& end, Args&&... args) {
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
    inline void evalExtFuncCallback(const PassivePosition& start, const PassivePosition& end, Args&&... args) {
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

    inline void evaluateInt(const Position& start, const Position& end) {
      primalValueCopy = (Real*)malloc(sizeof(Real) * primalsSize);
      memcpy(primalValueCopy, primals, sizeof(Real) * primalsSize);

      evaluateExtFunc(start, end);

      free(primalValueCopy);
    }

  public:

    /**
     * @brief Register a variable as an active variable.
     *
     * The index of the variable is set to a non zero number.
     * @param[inout] value The value which will be marked as an active variable.
     */
    inline void registerInput(ActiveReal<PrimalValueIndexTape<TapeTypes> >& value) {
      if(isActive()) {
        pushInputHandle(value.getValue(), value.getGradientData());
      }
    }

    /**
     * @brief Not needed in this implementation.
     *
     * @param[in] value Not used.
     */
    inline void registerOutput(ActiveReal<PrimalValueIndexTape<TapeTypes> >& value) {
      if(isActive() && value.getGradientData() != 0) {
        if(!IndexHandler::AssignNeedsStatement) {
          IndexType rhsIndex = value.getGradientData();

          indexHandler.assignIndex(value.getGradientData());
          pushCopyHandle(value.getValue(), value.getGradientData(), rhsIndex);
        }
      }
    }

    template<typename Stream = std::ostream>
    void printStatistics(Stream& out = std::cout) const {

      const std::string hLine = "-------------------------------------\n";

      size_t nChunksIndex  = indexVector.getNumChunks();
      size_t totalIndex    = indexVector.getDataSize();
      size_t sizeIndexEntry = sizeof(IndexType);
      double memoryUsedIndex = (double)totalIndex*(double)sizeIndexEntry* BYTE_TO_MB;
      double memoryAllocIndex= (double)nChunksIndex*(double)indexVector.getChunkSize()*(double)sizeIndexEntry* BYTE_TO_MB;

      size_t nChunksStmt  = stmtVector.getNumChunks();
      size_t totalStmt    = stmtVector.getDataSize();
      size_t sizeStmtEntry = sizeof(IndexType) + sizeof(Real) + sizeof(const ExpressionHandle<Real*, Real, IndexType>*) + sizeof(StatementInt);
      double memoryUsedStmt = (double)totalStmt*(double)sizeStmtEntry* BYTE_TO_MB;
      double memoryAllocStmt= (double)nChunksStmt*(double)stmtVector.getChunkSize()*(double)sizeStmtEntry* BYTE_TO_MB;

      size_t nChunksPassive  = constantValueVector.getNumChunks();
      size_t totalPassive    = constantValueVector.getDataSize();
      size_t sizePassiveEntry = sizeof(PassiveReal);
      double memoryUsedPassive = (double)totalPassive*(double)sizePassiveEntry* BYTE_TO_MB;
      double memoryAllocPassive= (double)nChunksPassive*(double)constantValueVector.getChunkSize()*(double)sizePassiveEntry* BYTE_TO_MB;

      size_t totalPrimal   = primalsSize;
      size_t sizePrimalEntry = sizeof(Real);
      double memoryAllocPrimal = (double)totalPrimal*(double)sizePrimalEntry* BYTE_TO_MB;

      out << hLine
          << "CoDi Tape Statistics (" << TapeTypes::tapeName << ")\n";
      printTapeBaseStatistics(out, hLine);
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

      printExtFuncStatistics(out, hLine);
    }

  private:
    inline void pushInputHandle(const Real& value, IndexType& index) {
      stmtVector.reserveItems(1);
      indexHandler.assignIndex(index);
      checkPrimalsSize();
      stmtVector.setDataAndMove(index, primals[index], &InputHandle, 0);

      primals[index] = value;

    }

    inline void pushCopyHandle(const Real& lhsValue, IndexType& lhsIndex, const IndexType& rhsIndex) {
      indexVector.reserveItems(1);
      indexVector.setDataAndMove(rhsIndex);

      stmtVector.reserveItems(1);
      checkPrimalsSize();
      stmtVector.setDataAndMove(lhsIndex, primals[lhsIndex], &CopyHandle, 0);

      primals[lhsIndex] = lhsValue;

    }

  };

  template <typename TapeTypes>
  const ExpressionHandle<typename TapeTypes::RealType*, typename TapeTypes::RealType, typename TapeTypes::IndexType> PrimalValueIndexTape<TapeTypes>::InputHandle(&PrimalValueIndexTape<TapeTypes>::inputHandleFunc, 0, 0);
  template <typename TapeTypes>
  const ExpressionHandle<typename TapeTypes::RealType*, typename TapeTypes::RealType, typename TapeTypes::IndexType> PrimalValueIndexTape<TapeTypes>::CopyHandle(&PrimalValueIndexTape<TapeTypes>::copyHandleFunc, 1, 0);

#define CREATE_PREACC_HANDLE(size) ExpressionHandle<typename TapeTypes::RealType*, typename TapeTypes::RealType, typename TapeTypes::IndexType>(&PrimalValueIndexTape<TapeTypes>::preaccHandleFunc<size>, size, size)

  template <typename TapeTypes>
  const ExpressionHandle<typename TapeTypes::RealType*, typename TapeTypes::RealType, typename TapeTypes::IndexType> PrimalValueIndexTape<TapeTypes>::PreaccHandles[MaxStatementIntSize] = {
    CREATE_PREACC_HANDLE(0),
    CREATE_PREACC_HANDLE(1),
    CREATE_PREACC_HANDLE(2),
    CREATE_PREACC_HANDLE(3),
    CREATE_PREACC_HANDLE(4),
    CREATE_PREACC_HANDLE(5),
    CREATE_PREACC_HANDLE(6),
    CREATE_PREACC_HANDLE(7),
    CREATE_PREACC_HANDLE(8),
    CREATE_PREACC_HANDLE(9),
    CREATE_PREACC_HANDLE(10),
    CREATE_PREACC_HANDLE(11),
    CREATE_PREACC_HANDLE(12),
    CREATE_PREACC_HANDLE(13),
    CREATE_PREACC_HANDLE(14),
    CREATE_PREACC_HANDLE(15),
    CREATE_PREACC_HANDLE(16),
    CREATE_PREACC_HANDLE(17),
    CREATE_PREACC_HANDLE(18),
    CREATE_PREACC_HANDLE(19),
    CREATE_PREACC_HANDLE(20),
    CREATE_PREACC_HANDLE(21),
    CREATE_PREACC_HANDLE(22),
    CREATE_PREACC_HANDLE(23),
    CREATE_PREACC_HANDLE(24),
    CREATE_PREACC_HANDLE(25),
    CREATE_PREACC_HANDLE(26),
    CREATE_PREACC_HANDLE(27),
    CREATE_PREACC_HANDLE(28),
    CREATE_PREACC_HANDLE(29),
    CREATE_PREACC_HANDLE(30),
    CREATE_PREACC_HANDLE(31),
    CREATE_PREACC_HANDLE(32),
    CREATE_PREACC_HANDLE(33),
    CREATE_PREACC_HANDLE(34),
    CREATE_PREACC_HANDLE(35),
    CREATE_PREACC_HANDLE(36),
    CREATE_PREACC_HANDLE(37),
    CREATE_PREACC_HANDLE(38),
    CREATE_PREACC_HANDLE(39),
    CREATE_PREACC_HANDLE(40),
    CREATE_PREACC_HANDLE(41),
    CREATE_PREACC_HANDLE(42),
    CREATE_PREACC_HANDLE(43),
    CREATE_PREACC_HANDLE(44),
    CREATE_PREACC_HANDLE(45),
    CREATE_PREACC_HANDLE(46),
    CREATE_PREACC_HANDLE(47),
    CREATE_PREACC_HANDLE(48),
    CREATE_PREACC_HANDLE(49),
    CREATE_PREACC_HANDLE(50),
    CREATE_PREACC_HANDLE(51),
    CREATE_PREACC_HANDLE(52),
    CREATE_PREACC_HANDLE(53),
    CREATE_PREACC_HANDLE(54),
    CREATE_PREACC_HANDLE(55),
    CREATE_PREACC_HANDLE(56),
    CREATE_PREACC_HANDLE(57),
    CREATE_PREACC_HANDLE(58),
    CREATE_PREACC_HANDLE(59),
    CREATE_PREACC_HANDLE(60),
    CREATE_PREACC_HANDLE(61),
    CREATE_PREACC_HANDLE(62),
    CREATE_PREACC_HANDLE(63),
    CREATE_PREACC_HANDLE(64),
    CREATE_PREACC_HANDLE(65),
    CREATE_PREACC_HANDLE(66),
    CREATE_PREACC_HANDLE(67),
    CREATE_PREACC_HANDLE(68),
    CREATE_PREACC_HANDLE(69),
    CREATE_PREACC_HANDLE(70),
    CREATE_PREACC_HANDLE(71),
    CREATE_PREACC_HANDLE(72),
    CREATE_PREACC_HANDLE(73),
    CREATE_PREACC_HANDLE(74),
    CREATE_PREACC_HANDLE(75),
    CREATE_PREACC_HANDLE(76),
    CREATE_PREACC_HANDLE(77),
    CREATE_PREACC_HANDLE(78),
    CREATE_PREACC_HANDLE(79),
    CREATE_PREACC_HANDLE(80),
    CREATE_PREACC_HANDLE(81),
    CREATE_PREACC_HANDLE(82),
    CREATE_PREACC_HANDLE(83),
    CREATE_PREACC_HANDLE(84),
    CREATE_PREACC_HANDLE(85),
    CREATE_PREACC_HANDLE(86),
    CREATE_PREACC_HANDLE(87),
    CREATE_PREACC_HANDLE(88),
    CREATE_PREACC_HANDLE(89),
    CREATE_PREACC_HANDLE(90),
    CREATE_PREACC_HANDLE(91),
    CREATE_PREACC_HANDLE(92),
    CREATE_PREACC_HANDLE(93),
    CREATE_PREACC_HANDLE(94),
    CREATE_PREACC_HANDLE(95),
    CREATE_PREACC_HANDLE(96),
    CREATE_PREACC_HANDLE(97),
    CREATE_PREACC_HANDLE(98),
    CREATE_PREACC_HANDLE(99),
    CREATE_PREACC_HANDLE(100),
    CREATE_PREACC_HANDLE(101),
    CREATE_PREACC_HANDLE(102),
    CREATE_PREACC_HANDLE(103),
    CREATE_PREACC_HANDLE(104),
    CREATE_PREACC_HANDLE(105),
    CREATE_PREACC_HANDLE(106),
    CREATE_PREACC_HANDLE(107),
    CREATE_PREACC_HANDLE(108),
    CREATE_PREACC_HANDLE(109),
    CREATE_PREACC_HANDLE(110),
    CREATE_PREACC_HANDLE(111),
    CREATE_PREACC_HANDLE(112),
    CREATE_PREACC_HANDLE(113),
    CREATE_PREACC_HANDLE(114),
    CREATE_PREACC_HANDLE(115),
    CREATE_PREACC_HANDLE(116),
    CREATE_PREACC_HANDLE(117),
    CREATE_PREACC_HANDLE(118),
    CREATE_PREACC_HANDLE(119),
    CREATE_PREACC_HANDLE(120),
    CREATE_PREACC_HANDLE(121),
    CREATE_PREACC_HANDLE(122),
    CREATE_PREACC_HANDLE(123),
    CREATE_PREACC_HANDLE(124),
    CREATE_PREACC_HANDLE(125),
    CREATE_PREACC_HANDLE(126),
    CREATE_PREACC_HANDLE(127),
    CREATE_PREACC_HANDLE(128),
    CREATE_PREACC_HANDLE(129),
    CREATE_PREACC_HANDLE(130),
    CREATE_PREACC_HANDLE(131),
    CREATE_PREACC_HANDLE(132),
    CREATE_PREACC_HANDLE(133),
    CREATE_PREACC_HANDLE(134),
    CREATE_PREACC_HANDLE(135),
    CREATE_PREACC_HANDLE(136),
    CREATE_PREACC_HANDLE(137),
    CREATE_PREACC_HANDLE(138),
    CREATE_PREACC_HANDLE(139),
    CREATE_PREACC_HANDLE(140),
    CREATE_PREACC_HANDLE(141),
    CREATE_PREACC_HANDLE(142),
    CREATE_PREACC_HANDLE(143),
    CREATE_PREACC_HANDLE(144),
    CREATE_PREACC_HANDLE(145),
    CREATE_PREACC_HANDLE(146),
    CREATE_PREACC_HANDLE(147),
    CREATE_PREACC_HANDLE(148),
    CREATE_PREACC_HANDLE(149),
    CREATE_PREACC_HANDLE(150),
    CREATE_PREACC_HANDLE(151),
    CREATE_PREACC_HANDLE(152),
    CREATE_PREACC_HANDLE(153),
    CREATE_PREACC_HANDLE(154),
    CREATE_PREACC_HANDLE(155),
    CREATE_PREACC_HANDLE(156),
    CREATE_PREACC_HANDLE(157),
    CREATE_PREACC_HANDLE(158),
    CREATE_PREACC_HANDLE(159),
    CREATE_PREACC_HANDLE(160),
    CREATE_PREACC_HANDLE(161),
    CREATE_PREACC_HANDLE(162),
    CREATE_PREACC_HANDLE(163),
    CREATE_PREACC_HANDLE(164),
    CREATE_PREACC_HANDLE(165),
    CREATE_PREACC_HANDLE(166),
    CREATE_PREACC_HANDLE(167),
    CREATE_PREACC_HANDLE(168),
    CREATE_PREACC_HANDLE(169),
    CREATE_PREACC_HANDLE(170),
    CREATE_PREACC_HANDLE(171),
    CREATE_PREACC_HANDLE(172),
    CREATE_PREACC_HANDLE(173),
    CREATE_PREACC_HANDLE(174),
    CREATE_PREACC_HANDLE(175),
    CREATE_PREACC_HANDLE(176),
    CREATE_PREACC_HANDLE(177),
    CREATE_PREACC_HANDLE(178),
    CREATE_PREACC_HANDLE(179),
    CREATE_PREACC_HANDLE(180),
    CREATE_PREACC_HANDLE(181),
    CREATE_PREACC_HANDLE(182),
    CREATE_PREACC_HANDLE(183),
    CREATE_PREACC_HANDLE(184),
    CREATE_PREACC_HANDLE(185),
    CREATE_PREACC_HANDLE(186),
    CREATE_PREACC_HANDLE(187),
    CREATE_PREACC_HANDLE(188),
    CREATE_PREACC_HANDLE(189),
    CREATE_PREACC_HANDLE(190),
    CREATE_PREACC_HANDLE(191),
    CREATE_PREACC_HANDLE(192),
    CREATE_PREACC_HANDLE(193),
    CREATE_PREACC_HANDLE(194),
    CREATE_PREACC_HANDLE(195),
    CREATE_PREACC_HANDLE(196),
    CREATE_PREACC_HANDLE(197),
    CREATE_PREACC_HANDLE(198),
    CREATE_PREACC_HANDLE(199),
    CREATE_PREACC_HANDLE(200),
    CREATE_PREACC_HANDLE(201),
    CREATE_PREACC_HANDLE(202),
    CREATE_PREACC_HANDLE(203),
    CREATE_PREACC_HANDLE(204),
    CREATE_PREACC_HANDLE(205),
    CREATE_PREACC_HANDLE(206),
    CREATE_PREACC_HANDLE(207),
    CREATE_PREACC_HANDLE(208),
    CREATE_PREACC_HANDLE(209),
    CREATE_PREACC_HANDLE(210),
    CREATE_PREACC_HANDLE(211),
    CREATE_PREACC_HANDLE(212),
    CREATE_PREACC_HANDLE(213),
    CREATE_PREACC_HANDLE(214),
    CREATE_PREACC_HANDLE(215),
    CREATE_PREACC_HANDLE(216),
    CREATE_PREACC_HANDLE(217),
    CREATE_PREACC_HANDLE(218),
    CREATE_PREACC_HANDLE(219),
    CREATE_PREACC_HANDLE(220),
    CREATE_PREACC_HANDLE(221),
    CREATE_PREACC_HANDLE(222),
    CREATE_PREACC_HANDLE(223),
    CREATE_PREACC_HANDLE(224),
    CREATE_PREACC_HANDLE(225),
    CREATE_PREACC_HANDLE(226),
    CREATE_PREACC_HANDLE(227),
    CREATE_PREACC_HANDLE(228),
    CREATE_PREACC_HANDLE(229),
    CREATE_PREACC_HANDLE(230),
    CREATE_PREACC_HANDLE(231),
    CREATE_PREACC_HANDLE(232),
    CREATE_PREACC_HANDLE(233),
    CREATE_PREACC_HANDLE(234),
    CREATE_PREACC_HANDLE(235),
    CREATE_PREACC_HANDLE(236),
    CREATE_PREACC_HANDLE(237),
    CREATE_PREACC_HANDLE(238),
    CREATE_PREACC_HANDLE(239),
    CREATE_PREACC_HANDLE(240),
    CREATE_PREACC_HANDLE(241),
    CREATE_PREACC_HANDLE(242),
    CREATE_PREACC_HANDLE(243),
    CREATE_PREACC_HANDLE(244),
    CREATE_PREACC_HANDLE(245),
    CREATE_PREACC_HANDLE(246),
    CREATE_PREACC_HANDLE(247),
    CREATE_PREACC_HANDLE(248),
    CREATE_PREACC_HANDLE(249),
    CREATE_PREACC_HANDLE(250),
    CREATE_PREACC_HANDLE(251),
    CREATE_PREACC_HANDLE(252),
    CREATE_PREACC_HANDLE(253),
    CREATE_PREACC_HANDLE(254),
    CREATE_PREACC_HANDLE(255)
  };

#undef CREATE_PREACC_HANDLE
}
