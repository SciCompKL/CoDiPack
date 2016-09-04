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

    /** @brief The termination of the vector sequence. */
    EmptyChunkVector emptyVector;

    #define TAPE_NAME PrimalValueIndexTape

    #define POSITION_TYPE typename TapeTypes::Position
    #define INDEX_HANDLER_TYPE IndexHandler
    #define RESET_FUNCTION_NAME resetExtFunc
    #define EVALUATE_FUNCTION_NAME evaluateInt
    #include "modules/tapeBaseModule.tpp"

    #define CHILD_VECTOR_TYPE EmptyChunkVector
    #define STMT_VECTOR_TYPE typename TapeTypes::StatementVector
    #define INDEX_VECTOR_TYPE typename TapeTypes::IndexVector
    #define CONSTANT_VECTOR_TYPE typename TapeTypes::ConstantValueVector
    #include "modules/primalValueModule.tpp"

    #define CHILD_VECTOR_TYPE ConstantValueVector
    #define CHILD_VECTOR_NAME constantValueVector
    #define VECTOR_TYPE typename TapeTypes::ExternalFunctionVector
    #include "modules/externalFunctionsModule.tpp"

    // TAPE_NAME is undefined at the end of the file

    Real* primalValueCopy;

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
      /* defined in the primalValueModule */stmtVector(DefaultChunkSize, emptyVector),
      /* defined in the primalValueModule */indexVector(DefaultChunkSize, stmtVector),
      /* defined in the primalValueModule */constantValueVector(DefaultChunkSize, indexVector),
      /* defined in the primalValueModule */primals(NULL),
      /* defined in the primalValueModule */primalsSize(0),
      /* defined in the primalValueModule */primalsIncr(DefaultSmallChunkSize),
      /* defined in externalFunctionsModule */extFuncVector(1000, constantValueVector),
      primalValueCopy(NULL) {
        // create the indices for the passive data
        for(size_t i = 0; i < MaxStatementIntSize; ++i) {
          indexHandler.createIndex();
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
    CODI_INLINE void clearAdjoints(const Position& start, const Position& end){
      CODI_UNUSED(start);
      CODI_UNUSED(end);
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

    CODI_INLINE void pushStmtData(IndexType& lhsIndex, const Real& rhsValue, const Handle& handle, const StatementInt& passiveVariableNumber) {
      indexHandler.assignIndex(lhsIndex);
      stmtVector.reserveItems(1);
      checkPrimalsSize();
      stmtVector.setDataAndMove(lhsIndex, primals[lhsIndex], handle, passiveVariableNumber);

      primals[lhsIndex] = rhsValue;
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
    CODI_INLINE void store(Real& lhsValue, IndexType& lhsIndex, const ActiveReal<PrimalValueIndexTape<TapeTypes> >& rhs) {
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
     * @brief Get the current position of the tape.
     *
     * The position can be used to reset the tape to that position or to
     * evaluate only parts of the tape.
     * @return The current position of the tape.
     */
    CODI_INLINE Position getPosition() const {
      return getExtFuncPosition();
    }

    /**
     * @brief Get the current position of the tape.
     *
     * The position can be used to reset the tape to that position or to
     * evaluate only parts of the tape.
     * @return The current position of the tape.
     */
    CODI_INLINE Position getZeroPosition() const {
      return getExtFuncZeroPosition();
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
    CODI_INLINE void evaluateStack(size_t& stmtPos, const size_t& endStmtPos, IndexType* lhsIndices, Real* storedPrimals, Handle* &statements, StatementInt* &passiveActiveReal, size_t& indexPos, IndexType* &indices, size_t& constantPos, PassiveReal* &constants, Real* primalVector) {
      while(stmtPos > endStmtPos) {
        --stmtPos;
        const IndexType& lhsIndex = lhsIndices[stmtPos];

        primalVector[lhsIndex] = storedPrimals[stmtPos];
        const GradientValue adj = adjoints[lhsIndex];
        adjoints[lhsIndex] = GradientValue();

        evaluateHandle(adj, statements[stmtPos], passiveActiveReal[stmtPos], indexPos, indices, constantPos, constants);
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
    CODI_INLINE void evalStmt(const StmtPosition& start, const StmtPosition& end, Args&&... args) {
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

    CODI_INLINE void evaluateInt(const Position& start, const Position& end) {
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
    CODI_INLINE void registerInput(ActiveReal<PrimalValueIndexTape<TapeTypes> >& value) {
      if(isActive()) {
        pushStmtData(value.getGradientData(), value.getValue(), &InputHandle, StatementInt(0));
      }
    }

    /**
     * @brief Not needed in this implementation.
     *
     * @param[in] value Not used.
     */
    CODI_INLINE void registerOutput(ActiveReal<PrimalValueIndexTape<TapeTypes> >& value) {
      if(isActive() && value.getGradientData() != 0) {
        if(!IndexHandler::AssignNeedsStatement) {
          IndexType rhsIndex = value.getGradientData();

          indexHandler.assignIndex(value.getGradientData());
          pushCopyHandle(value.getValue(), value.getGradientData(), rhsIndex);
        }
      }
    }

    /**
     * @brief Prints statistics about the tape on the screen or into a stream
     *
     * Prints information such as stored statements/adjoints and memory usage on screen or into
     * the stream when an argument is provided.
     *
     * @param[in,out] out  The information is written to the stream.
     *
     * @tparam Stream The type of the stream.
     */
    template<typename Stream = std::ostream>
    void printStatistics(Stream& out = std::cout) const {

      const std::string hLine = "-------------------------------------\n";

      out << hLine
          << "CoDi Tape Statistics (" << TapeTypes::tapeName << ")\n";
      printTapeBaseStatistics(out, hLine);
      printPrimalValueStatistics(out, hLine);
      printExtFuncStatistics(out, hLine);

    }
  };

  #include "modules/primalValueStaticModule.tpp"
  #undef TAPE_NAME

}
