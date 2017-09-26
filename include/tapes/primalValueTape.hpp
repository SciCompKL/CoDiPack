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


#pragma once

#include <cstddef>
#include <tuple>

#include "../activeReal.hpp"
#include "../expressionHandle.hpp"
#include "chunkVector.hpp"
#include "indices/linearIndexHandler.hpp"
#include "handles/functionHandleFactory.hpp"
#include "primalTapeExpressions.hpp"
#include "reverseTapeInterface.hpp"
#include "singleChunkVector.hpp"

namespace codi {

  /**
   * @brief Vector definition for the ChunkPrimalValueTape.
   *
   * The structure defines all vectors as chunk vectors.
   *
   * See PrimalValueTape for details.
   *
   * @tparam          Real  The type for the primal values.
   * @tparam  IndexHandler  The index handler for the managing of the indices. It has to be a index handler that assumes index reuse.
   * @tparam GradientValue  The type for the adjoint values. (Default: Same as the primal value.)
   * @tparam HandleFactory  The factory for the reverse interpretation of the expressions. Needs to implement the HandleFactoryInterface class.
   */
  template <typename Real, typename IndexHandler, typename GradientValue = Real, typename HandleFactory = FunctionHandleFactory<Real, typename IndexHandler::IndexType, Real> >
  struct ChunkPrimalValueTapeTypes {
    /** @brief The type for the primal values. */
    typedef Real RealType;
    /** @brief The handler for the indices. */
    typedef IndexHandler IndexHandlerType;
    /** @brief The type for the adjoint values. */
    typedef GradientValue GradientValueType;
    /** @brief The type for indices. */
    typedef typename IndexHandler::IndexType IndexType;

    /** @brief The type for the handle factory. */
    typedef HandleFactory HandleFactoryType;
    /** @brief The type for expression handles in the reverse evaluation. */
    typedef typename HandleFactory::Handle HandleType;

    /** @brief The data for each statement. */
    typedef Chunk2<HandleType, StatementInt> StatementChunk;
    /** @brief The chunk vector for the statement data. */
    typedef ChunkVector<StatementChunk, IndexHandler> StatementVector;

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

    /** @brief The name of the tape as a string. */
    constexpr static const char* tapeName = "ChunkPrimalValueTape";

  };

  /**
   * @brief Vector definition for the SimplePrimalValueTape.
   *
   * The structure defines all vectors as single chunk vectors.
   *
   * See PrimalValueTape for details.
   *
   * @tparam          Real  The type for the primal values.
   * @tparam  IndexHandler  The index handler for the managing of the indices. It has to be a index handler that assumes index reuse.
   * @tparam GradientValue  The type for the adjoint values. (Default: Same as the primal value.)
   * @tparam HandleFactory  The factory for the reverse interpretation of the expressions. Needs to implement the HandleFactoryInterface class.
   */
  template <typename Real, typename IndexHandler, typename GradientValue = Real, typename HandleFactory = FunctionHandleFactory<Real, typename IndexHandler::IndexType, Real> >
  struct SimplePrimalValueTapeTypes {
    /** @brief The type for the primal values. */
    typedef Real RealType;
    /** @brief The handler for the indices. */
    typedef IndexHandler IndexHandlerType;
    /** @brief The type for the adjoint values. */
    typedef GradientValue GradientValueType;
    /** @brief The type for indices. */
    typedef typename IndexHandler::IndexType IndexType;

    /** @brief The type for the handle factory. */
    typedef HandleFactory HandleFactoryType;
    /** @brief The type for expression handles in the reverse evaluation. */
    typedef typename HandleFactory::Handle HandleType;

    /** @brief The data for each statement. */
    typedef Chunk2<HandleType, StatementInt> StatementChunk;
    /** @brief The chunk vector for the statement data. */
    typedef SingleChunkVector<StatementChunk, IndexHandler> StatementVector;

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

    /** @brief The name of the tape as a string. */
    constexpr static const char* tapeName = "SimplePrimalValueTape";

  };

  /**
   * @brief A reverse AD tape that stores primal values for the reverse evaluation.
   *
   * The PrimalValueTape implements a fully featured ReverseTapeInterface. Depending on
   * the specified TapeTypes, new memory is automatically allocated or needs to be specified in advance.
   *
   * The current implementation uses 4 nested vectors
   * and the linear index handler as the terminator. The relation is
   *
   * externalFunctions -> constantValues -> indexData -> statements -> indexHandler
   *
   * The size of the tape can be set with the resize function,
   * the tape will allocate enough chunks such that the given data requirements will fit into the chunks.
   *
   * @tparam TapeTypes  All the types for the tape. Including the calculation type and the vector types.
   */
  template <typename TapeTypes>
  class PrimalValueTape : public ReverseTapeInterface<typename TapeTypes::RealType, typename TapeTypes::IndexHandlerType::IndexType, typename TapeTypes::GradientValueType, PrimalValueTape<TapeTypes>, typename TapeTypes::Position > {
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

    /** @brief The corresponding passive value to the real */
    typedef typename TypeTraits<Real>::PassiveReal PassiveReal;

    /** @brief The type of the handle factory */
    typedef typename TapeTypes::HandleFactoryType HandleFactory;

    /** @brief The type for expression handles in the reverse evaluation. */
    typedef typename TapeTypes::HandleType Handle;

    /** @brief The index handler for the active real's. */
    IndexHandler indexHandler;

    /** @brief Disables code path in CoDiPack that are optimized for Jacobi taping */
    static const bool AllowJacobiOptimization = false;

    #define TAPE_NAME PrimalValueTape

    #define POSITION_TYPE typename TapeTypes::Position
    #define INDEX_HANDLER_NAME indexHandler
    #define RESET_FUNCTION_NAME resetExtFunc
    #define EVALUATE_FUNCTION_NAME evaluateExtFunc
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

    #define ROOT_VECTOR extFuncVector
    #include "modules/ioModule.tpp"

    // TAPE_NAME is undefined at the end of the file


  public:
    /**
     * @brief Creates a tape with the size of zero for the data, statements and external functions.
     */
    PrimalValueTape() :
      indexHandler(MaxStatementIntSize - 1),
      /* defined in tapeBaseModule */adjoints(NULL),
      /* defined in tapeBaseModule */adjointsSize(0),
      /* defined in tapeBaseModule */active(false),
      /* defined in the primalValueModule */stmtVector(DefaultChunkSize, &indexHandler),
      /* defined in the primalValueModule */indexVector(DefaultChunkSize, &stmtVector),
      /* defined in the primalValueModule */constantValueVector(DefaultChunkSize, &indexVector),
      /* defined in the primalValueModule */primals(NULL),
      /* defined in the primalValueModule */primalsSize(0),
      /* defined in the primalValueModule */primalsIncr(DefaultSmallChunkSize),
      /* defined in externalFunctionsModule */extFuncVector(1000, &constantValueVector) {}

    /** @brief Tear down the tape. Delete all values from the modules */
    ~PrimalValueTape() {
      cleanTapeBase();
    }

    /**
     * @brief Swap the tape with an other tape.
     *
     * All data is exchanged between the tapes. The method performs the operation:
     *
     * T t = *this;
     * *this = other;
     * other = t;
     *
     */
    void swap(PrimalValueTape& other) {
      swapTapeBaseModule(other);
      swapPrimalValueModule(other);

      extFuncVector.swap(other.extFuncVector);
    }

    /**
     * @brief Sets all adjoint/gradients to zero.
     *
     * It has to hold start >= end.
     *
     * @param[in] start  The starting position for the reset of the vector.
     * @param[in]   end  The ending position for the reset of the vector.
     */
    CODI_INLINE void clearAdjoints(const Position& start, const Position& end) {

      IndexType startPos = min(end.inner.inner.inner.inner, adjointsSize - 1);
      IndexType endPos = min(start.inner.inner.inner.inner, adjointsSize - 1);

      for(IndexType i = startPos + 1; i <= endPos; ++i) {
        adjoints[i] = GradientValue();
      }
    }

    /**
     * @brief Set the size of the index and statement data and the primal vector.
     * @param[in] dataSize  The new size of the index vector.
     * @param[in] stmtSize  The new size of the statement vector.
     */
    void resize(const size_t& dataSize, const size_t& stmtSize) {
      indexVector.resize(dataSize);
      stmtVector.resize(stmtSize);

      resizePrimals(stmtSize + 1);
    }

    /**
     * @brief Pushes the handle to the statement vector and assigns a new index.
     *
     * The method also updates the value in the primal value vector.
     *
     * @param[in,out]          lhsIndex  The index of the lhs value. Will be renewed.
     * @param[in]              rhsValue  The value of the rhs. Is set in the primal value vector.
     * @param[in]                handle  The handle for the rhs expression.
     * @param[in] passiveVariableNumber  The number of passive values in the rhs.
     */
    CODI_INLINE void pushStmtData(IndexType& lhsIndex, const Real& rhsValue, const Handle& handle, const StatementInt& passiveVariableNumber) {
      stmtVector.reserveItems(1);
      stmtVector.setDataAndMove(handle, passiveVariableNumber);
      indexHandler.assignIndex(lhsIndex);

      checkPrimalsSize();
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
    CODI_INLINE void store(Real& lhsValue, IndexType& lhsIndex, const ActiveReal<PrimalValueTape<TapeTypes> >& rhs) {
      ENABLE_CHECK(OptTapeActivity, active){
        lhsIndex = rhs.getGradientData();
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
     * @param[in]       startAdjPos  The starting position for the adjoint evaluation.
     * @param[in]         endAdjPos  The ending position for the adjoint evaluation.
     * @param[in,out]       stmtPos  The current position in the statement data. It will decremented in the method.
     * @param[in]        statements  The vector with the handles for each statement.
     * @param[in] passiveActiveReal  The number passive values for each statement.
     * @param[in,out]      indexPos  The current position for the index data. It will decremented in the method.
     * @param[in]           indices  The indices for the arguments of the rhs.
     * @param[in,out]   constantPos  The current position in the constant data vector. It will decremented in the method.
     * @param[in]         constants  The constant values in the rhs expressions.
     */
    CODI_INLINE void evaluateStack(const size_t& startAdjPos, const size_t& endAdjPos, size_t& stmtPos, Handle* &statements, StatementInt* &passiveActiveReal, size_t& indexPos, IndexType* &indices, size_t& constantPos, PassiveReal* &constants) {
      size_t adjPos = startAdjPos;

      while(adjPos > endAdjPos) {
        const GradientValue& adj = adjoints[adjPos];
        --adjPos;
        --stmtPos;

        HandleFactory::template callHandle<PrimalValueTape<TapeTypes> >(statements[stmtPos], adj, passiveActiveReal[stmtPos], indexPos, indices, constantPos, constants, primals, adjoints);
      }
    }

    /**
     * @brief Evaluate a part of the statement vector.
     *
     * It has to hold start >= end.
     *
     * The function calls the evaluation method for the stack.
     *
     * @param[in] start  The starting point for the statement vector.
     * @param[in]   end  The ending point for the statement vector.
     * @param[in]  args  The arguments from the other vectors.
     *
     * @tparam Args  The types of the other arguments.
     */
    template<typename ... Args>
    CODI_INLINE void evalStmt(const StmtPosition& start, const StmtPosition& end, Args&&... args) {
      Handle* data;
      StatementInt* passiveActiveReal;
      size_t dataPos = start.data;
      auto curInnerPos = start.inner;
      for(size_t curChunk = start.chunk; curChunk > end.chunk; --curChunk) {
        stmtVector.getDataAtPosition(curChunk, 0, data, passiveActiveReal);

        auto endInnerPos = stmtVector.getInnerPosition(curChunk);
        evaluateStack(curInnerPos, endInnerPos, dataPos, data, passiveActiveReal, std::forward<Args>(args)...);

        codiAssert(dataPos == 0); // after a full chunk is evaluated, the data position needs to be zero

        curInnerPos = endInnerPos;

        dataPos = stmtVector.getChunkUsedData(curChunk - 1);
      }

      // Iterate over the reminder also covers the case if the start chunk and end chunk are the same
      stmtVector.getDataAtPosition(end.chunk, 0, data, passiveActiveReal);
      evaluateStack(curInnerPos, end.inner, dataPos, data, passiveActiveReal, std::forward<Args>(args)...);
    }

    /**
     * @brief Evaluate a part of the statement vector.
     *
     * It has to hold start >= end.
     *
     * The function calls the evaluation method for the constant value vector.
     *
     * @param[in] start  The starting point for the constant value vector.
     * @param[in]   end  The ending point for the constant value vector.
     * @param[in]  args  The arguments from the other vectors.
     *
     * @tparam Args  The types of the other arguments.
     */
    template<typename ... Args>
    CODI_INLINE void evalExtFuncCallback(const ConstantValuePosition& start, const ConstantValuePosition& end, Args&&... args) {
      evaluateConstantValues(start, end, std::forward<Args>(args)...);
    }


  public:

    /**
     * @brief Register a variable as an active variable.
     *
     * The index of the variable is set to a non zero number.
     * @param[in,out] value The value which will be marked as an active variable.
     */
    CODI_INLINE void registerInput(ActiveReal<PrimalValueTape<TapeTypes> >& value) {
      if(isActive()) {
        pushStmtData(value.getGradientData(), value.getValue(), HandleFactory::template createHandle<InputExpr<Real>, PrimalValueTape<TapeTypes> >(), StatementInt(0));
      }
    }

    /*
     * @brief It is ensured that each output variable has a unique index.
     *
     * @param[in] value  The value will have an unique index that is used by no other variable.
     */
    CODI_INLINE void registerOutput(ActiveReal<PrimalValueTape<TapeTypes> >& value) {
      if(isActive() && value.getGradientData() != 0) {
        IndexType rhsIndex = value.getGradientData();

        pushCopyHandle(value.getValue(), value.getGradientData(), rhsIndex);
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
