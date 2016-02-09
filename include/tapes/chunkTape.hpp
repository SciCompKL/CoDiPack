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

#pragma once

#include <iostream>
#include <iomanip>
#include <cstddef>
#include <tuple>

#include "../activeReal.hpp"
#include "chunk.hpp"
#include "chunkVector.hpp"
#include "externalFunctions.hpp"
#include "singleChunkVector.hpp"
#include "reverseTapeInterface.hpp"

/**
 * @brief Global namespace for CoDiPack - Code Differentiation Package
 */
namespace codi {

  /**
   * @brief Helper struct to define the nested chunk vectors for the ChunkTape.
   *
   * See ChunkTape for details.
   */
  template <typename Real, typename IndexHandler>
  struct ChunkTapeTypes {
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

  };

  /**
   * @brief Helper struct to define the nested vectors for the SimpleTape.
   *
   * See ChunkTape for details.
   */
  template <typename Real, typename IndexHandler>
  struct SimpleTapeTypes {
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

  };

  /**
   * @brief A tape which grows if more space is needed.
   *
   * The ChunkTape implements a fully featured ReverseTapeInterface in a most
   * user friendly fashion. The storage vectors of the tape are grown if the
   * tape runs out of space.
   *
   * This is handled by a nested definition of multiple ChunkVectors which hold
   * the different data vectors. The current implementation uses 3 ChunkVectors
   * and one terminator. The relation is
   *
   * externalFunctions -> statements -> jacobiData -> indexHandler.
   *
   * The size of the tape can be set with the resize function,
   * the tape will allocate enough chunks such that the given data requirements will fit into the chunks.
   *
   * @tparam        Real  The floating point type used in the ActiveReal.
   * @tparam   IndexType  The type for the indexing of the adjoint variables.
   * @tparam VectorTypes  The types for the data vectors.
   */
  template <typename Real, typename IndexHandler, typename VectorTypes>
  class ChunkTape : public ReverseTapeInterface<Real, typename IndexHandler::IndexType, ChunkTape<Real, IndexHandler, VectorTypes>, typename VectorTypes::Position > {
  public:

    #define TAPE_NAME ChunkTape

    typedef typename IndexHandler::IndexType IndexType;
    typedef IndexType GradientData;

    #define POSITION_TYPE typename VectorTypes::Position
    #define INDEX_HANDLER_TYPE IndexHandler
    #define RESET_FUNCTION_NAME resetExtFunc
    #define EVALUATE_FUNCTION_NAME evaluateExtFunc
    #include "modules/tapeBaseModule.tpp"

    #define CHILD_VECTOR_TYPE IndexHandler
    #define JACOBI_VECTOR_NAME jacobiVector
    #define VECTOR_TYPE typename VectorTypes::StatementVector
    #define STATEMENT_PUSH_FUNCTION_NAME pushStmtData
    #include "modules/statementModule.tpp"

    #define CHILD_VECTOR_TYPE StmtVector
    #define VECTOR_TYPE typename VectorTypes::JacobiVector
    #include "modules/jacobiModule.tpp"

    #define CHILD_VECTOR_TYPE JacobiVector
    #define CHILD_VECTOR_NAME jacobiVector
    #define VECTOR_TYPE typename VectorTypes::ExternalFunctionVector
    #include "modules/externalFunctionsModule.tpp"

    #undef TAPE_NAME

  public:
    /**
     * @brief Creates a tape with the default chunk sizes for the data, statements and
     * external functions defined in the configuration.
     */
    ChunkTape() :
      indexHandler(),
      adjoints(NULL),
      adjointsSize(0),
      active(false),
      stmtVector(DefaultChunkSize, indexHandler),
      jacobiVector(DefaultChunkSize, stmtVector),
      extFuncVector(1000, jacobiVector) {
    }

    inline void pushStmtData(const StatementInt& numberOfArguments, const IndexType& lhsIndex) {
      stmtVector.setDataAndMove(std::make_tuple(numberOfArguments));
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
    inline void store(Real& lhsValue, IndexType& lhsIndex, const ActiveReal<Real, ChunkTape<Real, IndexHandler, VectorTypes> >& rhs) {
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
    inline void clearAdjoints(const Position& start, const Position& end){
      for(IndexType i = end.inner.inner.inner; i <= start.inner.inner.inner; ++i) {
        adjoints[i] = 0.0;
      }
    }

    /**
     * @brief Get the current position of the tape.
     *
     * The position can be used to reset the tape to that position or to
     * evaluate only parts of the tape.
     * @return The current position of the tape.
     */
    inline Position getPosition() {
      return getExtFuncPosition();
    }


    /**
     * @brief Implementation of the AD stack evaluation.
     *
     * It has to hold startAdjPos >= endAdjPos.
     *
     * @param[in] startAdjPos The starting point in the expression evaluation.
     * @param[in]   endAdjPos The ending point in the expression evaluation.
     * @param[inout]  stmtPos The current position in the statement vector. This value is used in the next invocation of this method.
     * @param[in]  statements The pointer to the statement vector.
     * @param[inout]  dataPos The current position in the jacobi and index vector. This value is used in the next invocation of this method..
     * @param[in]    jacobies The pointer to the jacobi vector.
     * @param[in]     indices The pointer to the index vector
     */
    inline void evalStmtCallback(const size_t& startAdjPos, const size_t& endAdjPos, size_t& stmtPos, StatementInt* &statements, size_t& dataPos, Real* &jacobies, IndexType* &indices) {
      size_t adjPos = startAdjPos;

      while(adjPos > endAdjPos) {
        const Real& adj = adjoints[adjPos];
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
     * @param[in] start The starting point for the statement vector.
     * @param[in]   end The ending point for the statement vector.
     */
    template<typename ... Args>
    inline void evaluateStmt(const StmtPosition& start, const StmtPosition& end, Args&&... args) {
      StatementInt* statementData;
      size_t dataPos = start.data;
      StmtChildPosition curInnerPos = start.inner;
      for(size_t curChunk = start.chunk; curChunk > end.chunk; --curChunk) {
        std::tie(statementData) = stmtVector.getDataAtPosition(curChunk, 0);

        StmtChildPosition endInnerPos = stmtVector.getInnerPosition(curChunk);
        evalStmtCallback(curInnerPos, endInnerPos, dataPos, statementData, std::forward<Args>(args)...);

        curInnerPos = endInnerPos;

        dataPos = stmtVector.getChunkUsedData(curChunk - 1);
      }

      // Iterate over the reminder also covers the case if the start chunk and end chunk are the same
      std::tie(statementData) = stmtVector.getDataAtPosition(end.chunk, 0);
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
    inline void evalJacobiesCallback(const StmtPosition& start, const StmtPosition& end, size_t& dataPos, Real* &jacobies, IndexType* &indices) {
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
    inline void evalExtFuncCallback(const JacobiPosition& start, const JacobiPosition& end) {
      evaluateJacobies(start, end);
    }

  public:

    /**
     * @brief Register a variable as an active variable.
     *
     * The index of the variable is set to a non zero number.
     * @param[inout] value The value which will be marked as an active variable.
     */
    inline void registerInput(ActiveReal<Real, ChunkTape<Real, IndexHandler, VectorTypes> >& value) {
      stmtVector.reserveItems(1);
      stmtVector.setDataAndMove(std::make_tuple((StatementInt)0));

      value.getGradientData() = indexHandler.createIndex();
    }

    /**
     * @brief Register the value as an ouput value.
     *
     * The method ensures that each output value has an unique index. This is done
     * by performing the trivial operation value *= 1.0
     *
     * @param[in] value A new index is assigned.
     */
    inline void registerOutput(ActiveReal<Real, ChunkTape<Real, IndexHandler, VectorTypes> >& value) {
      value = 1.0 * value;
    }

    /**
     * @brief Prints statistics about the tape on the screen
     *
     * Prints information such as stored statements/adjoints and memory usage on screen.
     */
    void printStatistics(){

      std::cout << std::endl
                << "-------------------------------------" << std::endl
                << "CoDi Tape Statistics (ChunkTape)"      << std::endl
                << "-------------------------------------" << std::endl;
      printTapeBaseStatistics();
      printStmtStatistics();
      printJacobiStatistics();
      printExtFuncStatistics();
      indexHandler.printStatistics();
      std::cout << std::endl;

    }

  };
}
