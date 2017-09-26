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

#include <iostream>
#include <iomanip>
#include <cstddef>
#include <tuple>

#include "../activeReal.hpp"
#include "../typeFunctions.hpp"
#include "chunk.hpp"
#include "chunkVector.hpp"
#include "externalFunctions.hpp"
#include "reverseTapeInterface.hpp"
#include "singleChunkVector.hpp"

/**
 * @brief Global namespace for CoDiPack - Code Differentiation Package
 */
namespace codi {

  /**
   * @brief Vector definition for the ChunkTape.
   *
   * The structure defines all vectors as chunk vectors.
   *
   * See JacobiTape for details.
   *
   * @tparam Real  The type for the primal values.
   * @tparam IndexHandler  The index handler for the managing of the indices. It has to be a index handler that assumes index reuse.
   * @tparam GradientValue  The type for the adjoint values. (Default: Same as the primal value.)
   */
  template <typename Real, typename IndexHandler, typename GradientValue = Real>
  struct ChunkTapeTypes {
    /** @brief The type for the primal values. */
    typedef Real RealType;
    /** @brief The handler for the indices. */
    typedef IndexHandler IndexHandlerType;
    /** @brief The type for the adjoint values. */
    typedef GradientValue GradientValueType;

    /** @brief The data for each statement. */
    typedef Chunk1<StatementInt> StatementChunk;
    /** @brief The chunk vector for the statement data. */
    typedef ChunkVector<StatementChunk, IndexHandler> StatementVector;

    /** @brief The data for the jacobies of each statement */
    typedef Chunk2< Real, typename IndexHandler::IndexType> JacobiChunk;
    /** @brief The chunk vector for the jacobi data. */
    typedef ChunkVector<JacobiChunk, StatementVector> JacobiVector;

    /** @brief The data for the external functions. */
    typedef Chunk2<ExternalFunction,typename JacobiVector::Position> ExternalFunctionChunk;
    /** @brief The chunk vector for the external  function data. */
    typedef ChunkVector<ExternalFunctionChunk, JacobiVector> ExternalFunctionVector;

    /** @brief The position for all the different data vectors. */
    typedef typename ExternalFunctionVector::Position Position;

    /** @brief The name of the tape as a string. */
    constexpr static const char* tapeName = "ChunkTape";

  };

  /**
   * @brief Vector definition for the SimpleTape.
   *
   * The structure defines all vectors as single chunk vectors.
   *
   * See JacobiTape for details.
   *
   * @tparam Real  The type for the primal values.
   * @tparam IndexHandler  The index handler for the managing of the indices. It has to be a index handler that assumes index reuse.
   * @tparam GradientValue  The type for the adjoint values. (Default: Same as the primal value.)
   */
  template <typename Real, typename IndexHandler, typename GradientValue = Real>
  struct SimpleTapeTypes {
    /** @brief The type for the primal values. */
    typedef Real RealType;
    /** @brief The handler for the indices. */
    typedef IndexHandler IndexHandlerType;
    /** @brief The type for the adjoint values. */
    typedef GradientValue GradientValueType;

    /** @brief The data for each statement. */
    typedef Chunk1<StatementInt> StatementChunk;
    /** @brief The chunk vector for the statement data. */
    typedef SingleChunkVector<StatementChunk, IndexHandler> StatementVector;

    /** @brief The data for the jacobies of each statement */
    typedef Chunk2< Real, typename IndexHandler::IndexType> JacobiChunk;
    /** @brief The chunk vector for the jacobi data. */
    typedef SingleChunkVector<JacobiChunk, StatementVector> JacobiVector;

    /** @brief The data for the external functions. */
    typedef Chunk2<ExternalFunction,typename JacobiVector::Position> ExternalFunctionChunk;
    /** @brief The chunk vector for the external  function data. */
    typedef SingleChunkVector<ExternalFunctionChunk, JacobiVector> ExternalFunctionVector;

    /** @brief The position for all the different data vectors. */
    typedef typename ExternalFunctionVector::Position Position;

    /** @brief The name of the tape as a string. */
    constexpr static const char* tapeName = "SimpleTape";

  };

  /**
   * @brief A reverse AD tape that stores Jacobie values for the reverse evaluation.
   *
   * The JacobiTape implements a fully featured ReverseTapeInterface. Depending on
   * the specified TapeTypes new memory, is automatically allocated or needs to be specified in advance.
   *
   * The current implementation uses 3 nested vectors
   * and the linear index handler as the terminator. The relation is
   *
   * externalFunctions -> jacobiData -> statements -> indexHandler
   *
   * The size of the tape can be set with the resize function,
   * the tape will allocate enough chunks such that the given data requirements will fit into the chunks.
   *
   * @tparam TapeTypes  All the types for the tape. Including the calculation type and the vector types.
   */
  template <typename TapeTypes>
  class JacobiTape final : public ReverseTapeInterface<typename TapeTypes::RealType, typename TapeTypes::IndexHandlerType::IndexType, typename TapeTypes::GradientValueType, JacobiTape<TapeTypes>, typename TapeTypes::Position > {
  public:

    /** @brief The type for the primal values. */
    typedef typename TapeTypes::RealType Real;
    /** @brief The type for the adjoint values. */
    typedef typename TapeTypes::GradientValueType GradientValue;
    /** @brief The type for the index handler. */
    typedef typename TapeTypes::IndexHandlerType IndexHandler;

    /** @brief The type for the indices that are used for the identification of the adjoint variables. */
    typedef typename IndexHandler::IndexType IndexType;
    /** @brief The gradient data is just the index type. */
    typedef IndexType GradientData;

    /** @brief The corresponding passive value to the real */
    typedef typename TypeTraits<Real>::PassiveReal PassiveReal;

    /** @brief The index handler for the active real's. */
    IndexHandler indexHandler;

    /** @brief Enables code path in CoDiPack that are optimized for Jacobi taping */
    static const bool AllowJacobiOptimization = true;

    // The class name of the tape. Required by the modules.
    #define TAPE_NAME JacobiTape

    #define POSITION_TYPE typename TapeTypes::Position
    #define INDEX_HANDLER_NAME indexHandler
    #define RESET_FUNCTION_NAME resetExtFunc
    #define EVALUATE_FUNCTION_NAME evaluateExtFunc
    #include "modules/tapeBaseModule.tpp"

    #define CHILD_VECTOR_TYPE IndexHandler
    #define JACOBI_VECTOR_NAME jacobiVector
    #define VECTOR_TYPE typename TapeTypes::StatementVector
    #define STATEMENT_PUSH_FUNCTION_NAME pushStmtData
    #include "modules/statementModule.tpp"

    #define CHILD_VECTOR_TYPE StmtVector
    #define VECTOR_TYPE typename TapeTypes::JacobiVector
    #include "modules/jacobiModule.tpp"

    #define CHILD_VECTOR_TYPE JacobiVector
    #define CHILD_VECTOR_NAME jacobiVector
    #define VECTOR_TYPE typename TapeTypes::ExternalFunctionVector
    #include "modules/externalFunctionsModule.tpp"

    #define ROOT_VECTOR extFuncVector
    #include "modules/ioModule.tpp"

    #undef TAPE_NAME

  public:
    /**
     * @brief Creates a tape with the default chunk sizes for the data, statements and
     * external functions defined in the configuration.
     */
    JacobiTape() :
      indexHandler(0),
      /* defined in tapeBaseModule */adjoints(NULL),
      /* defined in tapeBaseModule */adjointsSize(0),
      /* defined in tapeBaseModule */active(false),
      /* defined in statementModule */stmtVector(DefaultChunkSize, &indexHandler),
      /* defined in jacobiModule */jacobiVector(DefaultChunkSize, &stmtVector),
      /* defined in externalFunctionsModule */extFuncVector(1000, &jacobiVector) {
    }

    /** @brief Tear down the tape. Delete all values from the modules */
    ~JacobiTape() {
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
    void swap(JacobiTape& other) {
      swapTapeBaseModule(other);

      extFuncVector.swap(other.extFuncVector);
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
    CODI_INLINE void store(Real& lhsValue, IndexType& lhsIndex, const ActiveReal<JacobiTape<TapeTypes> >& rhs) {
      ENABLE_CHECK (OptTapeActivity, active){
        lhsIndex = rhs.getGradientData();
      } else {
        indexHandler.freeIndex(lhsIndex);
      }
      lhsValue = rhs.getValue();
    }

    /**
     * @brief Set the size of the jacobi and statement data.
     *
     * The tape will allocate enough chunks such that the given data
     * sizes will fit into the chunk vectors.
     *
     * @param[in]      dataSize  The new size of the jacobi data.
     * @param[in] statementSize  The new size of the statement data.
     */
    void resize(const size_t& dataSize, const size_t& statementSize) {
      resizeJacobi(dataSize);
      resizeStmt(statementSize);
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
      IndexType startPos = min(end.inner.inner.inner, adjointsSize - 1);
      IndexType endPos = min(start.inner.inner.inner, adjointsSize - 1);

      for(IndexType i = startPos + 1; i <= endPos; ++i) {
        adjoints[i] = GradientValue();
      }
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
     * @brief Get the initial position of the tape.
     *
     * The position can be used to reset the tape to that position or to
     * evaluate only parts of the tape.
     * @return The initial position of the tape.
     */
    CODI_INLINE Position getZeroPosition() const {
      return getExtFuncZeroPosition();
    }


  private:

    /**
     * @brief The callback method for the push of statement data.
     *
     * The method is called by the statement module to push the
     * statements on the tape.
     *
     * @param[in] numberOfArguments  The number of arguments in the statements that have been pushed as jacobies.
     * @param[in]          lhsIndex  The index of the lhs value of the operation.
     */
    CODI_INLINE void pushStmtData(const StatementInt& numberOfArguments, const IndexType& lhsIndex) {
      CODI_UNUSED(lhsIndex);

      stmtVector.setDataAndMove(numberOfArguments);
    }

    /**
     * @brief Implementation of the AD stack evaluation.
     *
     * It has to hold startAdjPos >= endAdjPos.
     *
     * @param[in] startAdjPos The starting point in the expression evaluation.
     * @param[in]   endAdjPos The ending point in the expression evaluation.
     * @param[in,out] stmtPos The current position in the statement vector. This value is used in the next invocation of this method.
     * @param[in]  statements The pointer to the statement vector.
     * @param[in,out] dataPos The current position in the jacobi and index vector. This value is used in the next invocation of this method..
     * @param[in]    jacobies The pointer to the jacobi vector.
     * @param[in]     indices The pointer to the index vector
     */
    CODI_INLINE void evalStmtCallback(const size_t& startAdjPos, const size_t& endAdjPos, size_t& stmtPos, StatementInt* &statements, size_t& dataPos, Real* &jacobies, IndexType* &indices) {
      size_t adjPos = startAdjPos;

      while(adjPos > endAdjPos) {
        const GradientValue& adj = adjoints[adjPos];
        --adjPos;
        --stmtPos;

        incrementAdjoints(adj, adjoints, statements[stmtPos], dataPos, jacobies, indices);
      }
    }

    /**
     * @brief Evaluate a part of the statement vector.
     *
     * It has to hold start >= end.
     *
     * The function calls the evaluation method for the jacobi vector.
     *
     * @param[in] start  The starting point for the statement vector.
     * @param[in]   end  The ending point for the statement vector.
     * @param[in]  args  The arguments from the other vectors.
     *
     * @tparam Args  The types of the other arguments.
     */
    template<typename ... Args>
    CODI_INLINE void evaluateStmt(const StmtPosition& start, const StmtPosition& end, Args&&... args) {
      StatementInt* statementData;
      size_t dataPos = start.data;
      StmtChildPosition curInnerPos = start.inner;
      for(size_t curChunk = start.chunk; curChunk > end.chunk; --curChunk) {
        stmtVector.getDataAtPosition(curChunk, 0, statementData);

        StmtChildPosition endInnerPos = stmtVector.getInnerPosition(curChunk);
        evalStmtCallback(curInnerPos, endInnerPos, dataPos, statementData, std::forward<Args>(args)...);

        codiAssert(dataPos == 0); // after a full chunk is evaluated, the data position needs to be zero

        curInnerPos = endInnerPos;

        dataPos = stmtVector.getChunkUsedData(curChunk - 1);
      }

      // Iterate over the reminder also covers the case if the start chunk and end chunk are the same
      stmtVector.getDataAtPosition(end.chunk, 0, statementData);
      evalStmtCallback(curInnerPos, end.inner, dataPos, statementData, std::forward<Args>(args)...);
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
    CODI_INLINE void evalJacobiesCallback(const StmtPosition& start, const StmtPosition& end, size_t& dataPos, Real* &jacobies, IndexType* &indices) {
      evaluateStmt(start, end, dataPos, jacobies, indices);
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
    CODI_INLINE void evalExtFuncCallback(const JacobiPosition& start, const JacobiPosition& end) {
      evaluateJacobies(start, end);
    }

  public:

    /**
     * @brief Register a variable as an active variable.
     *
     * The index of the variable is set to a non zero number.
     * @param[in,out] value The value which will be marked as an active variable.
     */
    CODI_INLINE void registerInput(ActiveReal<JacobiTape<TapeTypes> >& value) {
      stmtVector.reserveItems(1);
      stmtVector.setDataAndMove((StatementInt)0);

      value.getGradientData() = indexHandler.createIndex();
    }

    /**
     * @brief Register the value as an output value.
     *
     * The method ensures that each output value has an unique index. This is done
     * by performing the trivial operation value *= 1.0
     *
     * @param[in] value A new index is assigned.
     */
    CODI_INLINE void registerOutput(ActiveReal<JacobiTape<TapeTypes> >& value) {
      value = 1.0 * value;
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
      printStmtStatistics(out, hLine);
      printJacobiStatistics(out, hLine);
      printExtFuncStatistics(out, hLine);
      indexHandler.printStatistics(out, hLine);

    }
  };
}
