/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2018 Chair for Scientific Computing (SciComp), TU Kaiserslautern
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
#include <cstddef>
#include <tuple>

#include "../activeReal.hpp"
#include "../typeFunctions.hpp"
#include "chunk.hpp"
#include "chunkVector.hpp"
#include "externalFunctions.hpp"
#include "reverseTapeInterface.hpp"
#include "singleChunkVector.hpp"
#include "../tapeTypes.hpp"
#include "../tools/tapeValues.hpp"

/**
 * @brief Global namespace for CoDiPack - Code Differentiation Package
 */
namespace codi {

  /**
   * @brief Vector definition for the ChunkIndexTape.
   *
   * The structure defines all vectors as chunk vectors.
   *
   * See JacobiIndexTape for details.
   *
   * @tparam        RTT  The basic type definitions for the tape. Need to define everything from ReverseTapeTypes.
   * @tparam DataVector  The data manager for the chunks. Needs to implement a ChunkVector interface.
   */
  template <typename RTT, template<typename, typename> class DataVector>
  struct JacobiIndexTapeTypes {

    CODI_INLINE_REVERSE_TAPE_TYPES(RTT)

    /** @brief The tape type structure, that defines the basic types. */
    typedef RTT BaseTypes;

    /** @brief The data for each statement. */
    typedef Chunk2<StatementInt, typename IndexHandler::Index> StatementChunk;
    /** @brief The chunk vector for the statement data. */
    typedef ChunkVector<StatementChunk, EmptyChunkVector> StatementVector;

    /** @brief The data for the jacobies of each statement */
    typedef Chunk2< Real, typename IndexHandler::Index> JacobiChunk;
    /** @brief The chunk vector for the jacobi data. */
    typedef ChunkVector<JacobiChunk, StatementVector> JacobiVector;

    /** @brief The data for the external functions. */
    typedef Chunk2<ExternalFunction,typename JacobiVector::Position> ExternalFunctionChunk;
    /** @brief The chunk vector for the external  function data. */
    typedef ChunkVector<ExternalFunctionChunk, JacobiVector> ExternalFunctionVector;

    /** @brief The position for all the different data vectors. */
    typedef typename ExternalFunctionVector::Position Position;

    /** @brief The name of the tape as a string. */
    constexpr static const char* tapeName = "JacobiIndexTape";

  };

  /**
   * @brief A reverse AD tape that stores Jacobie values for the reverse evaluation.
   *
   * The JacobiIndexTape implements a fully featured ReverseTapeInterface. Depending on
   * the specified TapeTypes, new memory is automatically allocated or needs to be specified in advance.
   *
   * The current implementation uses 3 nested vectors
   * and an empty terminator. The relation is
   *
   * externalFunctions -> jacobiData -> statements
   *
   * The size of the tape can be set with the resize function,
   * the tape will allocate enough chunks such that the given data requirements will fit into the chunks.
   *
   * The tape also uses the specified index manager to reuse the indices that are deleted.
   * That means that ActiveReal's which use this tape need to be copied by usual means and deleted after
   * they are no longer used. No c-like memory operations like memset and memcpy should be applied
   * to these types.
   *
   * @tparam TapeTypes  All the types for the tape. Including the calculation type and the vector types.
   */
  template <typename TapeTypes>
  class JacobiIndexTape final : public ReverseTapeInterface<typename TapeTypes::Real, typename TapeTypes::Index, typename TapeTypes::GradientValue, JacobiIndexTape<TapeTypes>, typename TapeTypes::Position > {
  public:

    CODI_INLINE_REVERSE_TAPE_TYPES(TapeTypes::BaseTypes)

    /** @brief The tape type structure, that defines the basic types. */
    typedef typename TapeTypes::BaseTypes BaseTypes;

    /** @brief The gradient data for this tape is the type of the indices. */
    typedef Index GradientData;

    /** @brief The counter for the current expression. */
    EmptyChunkVector emptyVector;

    /** @brief The index handler for the active real's. */
    static typename TapeTypes::IndexHandler indexHandler;

    /** @brief Enables code path in CoDiPack that are optimized for Jacobi taping */
    static const bool AllowJacobiOptimization = true;

    /** @brief This tape requires no primal value handling. */
    static const bool RequiresPrimalReset = false;

    // The class name of the tape. Required by the modules.
    #define TAPE_NAME JacobiIndexTape

    #define POSITION_TYPE typename TapeTypes::Position
    #define INDEX_HANDLER_NAME indexHandler
    #define RESET_FUNCTION_NAME resetExtFunc
    #define EVALUATE_FUNCTION_NAME evaluateInt
    #define EVALUATE_FORWARD_FUNCTION_NAME evaluateForwardInt
    #define EVALUATE_PRIMAL_FUNCTION_NAME evaluatePrimalStub
    #include "modules/tapeBaseModule.tpp"

    #define CHILD_VECTOR_TYPE EmptyChunkVector
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
    JacobiIndexTape() :
      emptyVector(),
      /* defined in tapeBaseModule */adjoints(NULL),
      /* defined in tapeBaseModule */adjointsSize(0),
      /* defined in tapeBaseModule */active(false),
      /* defined in statementModule */stmtVector(DefaultChunkSize, &emptyVector),
      /* defined in jacobiModule */jacobiVector(DefaultChunkSize, &stmtVector),
      /* defined in externalFunctionsModule */extFuncVector(1000, &jacobiVector) {
    }

    /** @brief Tear down the tape. Delete all values from the modules */
    ~JacobiIndexTape() {
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
    void swap(JacobiIndexTape& other) {
      swapTapeBaseModule(other);

      // the index handler is not swapped because the indices of the program state need to stay valid

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
    CODI_INLINE void store(Real& lhsValue, Index& lhsIndex, const ActiveReal<JacobiIndexTape<TapeTypes> >& rhs) {
      ENABLE_CHECK (OptTapeActivity, active){
        ENABLE_CHECK(OptCheckZeroIndex, 0 != rhs.getGradientData()) {
          indexHandler.copyIndex(lhsIndex, rhs.getGradientData());

          if(IndexHandler::AssignNeedsStatement) {
            stmtVector.reserveItems(1);
            jacobiVector.reserveItems(1);
            jacobiVector.setDataAndMove(PassiveReal(1.0), rhs.getGradientData());
            stmtVector.setDataAndMove((StatementInt)1, lhsIndex);
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
      if(NULL != adjoints) {
        auto clearFunc = [this] (StatementInt* stmtSize, Index* index) {
          CODI_UNUSED(stmtSize);

          if(*index < adjointsSize) {
            this->adjoints[*index] = GradientValue();
          }
        };

        stmtVector.forEachReverse(start.inner.inner, end.inner.inner, clearFunc);
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
    CODI_INLINE void pushStmtData(const StatementInt& numberOfArguments, const Index& lhsIndex) {
      stmtVector.setDataAndMove(numberOfArguments, lhsIndex);
    }

    /**
     * @brief Implementation of the AD stack evaluation.
     *
     * It has to hold startAdjPos >= endAdjPos.
     *
     * @param[in,out]      adjointData  The vector of the adjoint variables.
     * @param[in,out]          dataPos The current position in the jacobi and index vector. This value is used in the next invocation of this method..
     * @param[in]           endDataPos The end position in the jacobi and index vector.
     * @param[in]             jacobies The pointer to the jacobies of the rhs arguments.
     * @param[in]              indices The pointer the indices of the rhs arguments.
     * @param[in,out]          stmtPos The starting point in the expression evaluation. The index is decremented.
     * @param[in]           endStmtPos The ending point in the expression evaluation.
     * @param[in]    numberOfArguments The pointer to the number of arguments of the statement.
     * @param[in]           lhsIndices The pointer the indices of the lhs.
     *
     * @tparam AdjointData The data for the adjoint vector it needs to support add, multiply and comparison operations.
     */
    template<typename AdjointData>
    CODI_INLINE void evaluateStackReverse(AdjointData* adjointData,
                                          size_t& dataPos, const size_t& endDataPos, Real* &jacobies, Index* &indices,
                                          size_t& stmtPos, const size_t& endStmtPos, StatementInt* &numberOfArguments,
                                          Index* lhsIndices) {

      CODI_UNUSED(endDataPos);

      while(stmtPos > endStmtPos) {
        --stmtPos;
        const Index& lhsIndex = lhsIndices[stmtPos];
        const AdjointData adj = adjointData[lhsIndex];
        adjointData[lhsIndex] = GradientValue();

#if CODI_AdjointHandle_Jacobi_Reverse
        handleReverseEval(adj, lhsIndex);
#endif

        incrementAdjoints(adj, adjointData, numberOfArguments[stmtPos], dataPos, jacobies, indices);
      }
    }

    /**
     * @brief Evaluate the stack in reverse order.
     *
     * It has to hold start >= end.
     *
     * @param[in]           start  The start point for the evaluation.
     * @param[in]             end  The end point for the evaluation.
     * @param[in,out] adjointData  The vector of the adjoint variables.
     *
     * @tparam AdjointData The data for the adjoint vector it needs to support add, multiply and comparison operations.
     */
    template<typename AdjointData>
    CODI_INLINE void evaluateInt(const Position& start, const Position& end, AdjointData* adjointData) {

      auto evalFunc = [this] (AdjointData* adjointData,
          size_t& dataPos, const size_t& endDataPos, Real* &jacobies, Index* &indices,
          size_t& stmtPos, const size_t& endStmtPos, StatementInt* &statements, Index* lhsIndices) {
        evaluateStackReverse<AdjointData>(adjointData, dataPos, endDataPos, jacobies, indices,
                                      stmtPos, endStmtPos, statements, lhsIndices);
      };
      auto reverseFunc = &JacobiVector::template evaluateReverse<decltype(evalFunc), AdjointData*&>;

      AdjointInterfaceImpl<Real, Index, AdjointData> interface(adjointData);

      evaluateExtFunc(start, end, reverseFunc, jacobiVector, &interface, evalFunc, adjointData);
    }

    /**
     * @brief Implementation of the AD stack evaluation.
     *
     * It has to hold startAdjPos <= endAdjPos.
     *
     * @param[in,out]      adjointData  The vector of the adjoint variables.
     * @param[in,out]          dataPos The current position in the jacobi and index vector. This value is used in the next invocation of this method..
     * @param[in]           endDataPos The end position in the jacobi and index vector.
     * @param[in]             jacobies The pointer to the jacobies of the rhs arguments.
     * @param[in]              indices The pointer the indices of the rhs arguments.
     * @param[in,out]          stmtPos The starting point in the expression evaluation. The index is decremented.
     * @param[in]           endStmtPos The ending point in the expression evaluation.
     * @param[in]    numberOfArguments The pointer to the number of arguments of the statement.
     * @param[in]           lhsIndices The pointer the indices of the lhs.
     *
     * @tparam AdjointData The data for the adjoint vector it needs to support add, multiply and comparison operations.
     */
    template<typename AdjointData>
    CODI_INLINE void evaluateStackForward(AdjointData* adjointData,
                                          size_t& dataPos, const size_t& endDataPos, Real* &jacobies, Index* &indices,
                                          size_t& stmtPos, const size_t& endStmtPos, StatementInt* &numberOfArguments,
                                          Index* lhsIndices) {
      CODI_UNUSED(endDataPos);

      while(stmtPos < endStmtPos) {
        const Index& lhsIndex = lhsIndices[stmtPos];
        AdjointData adj = AdjointData();

        incrementTangents(adj, adjointData, numberOfArguments[stmtPos], dataPos, jacobies, indices);
        adjointData[lhsIndex] = adj;

        ++stmtPos;
      }
    }

    /**
     * @brief Evaluate the stack in forward order.
     *
     * It has to hold start <= end.
     *
     * @param[in]           start  The start point for the evaluation.
     * @param[in]             end  The end point for the evaluation.
     * @param[in,out] adjointData  The vector of the adjoint variables.
     *
     * @tparam AdjointData The data for the adjoint vector it needs to support add, multiply and comparison operations.
     */
    template<typename AdjointData>
    CODI_INLINE void evaluateForwardInt(const Position& start, const Position& end, AdjointData* adjointData) {

      auto evalFunc = [this] (AdjointData* adjointData,
          size_t& dataPos, const size_t& endDataPos, Real* &jacobies, Index* &indices,
          size_t& stmtPos, const size_t& endStmtPos, StatementInt* &statements, Index* lhsIndices) {
        evaluateStackForward<AdjointData>(adjointData, dataPos, endDataPos, jacobies, indices,
                                      stmtPos, endStmtPos, statements, lhsIndices);
      };
      auto forwardFunc = &JacobiVector::template evaluateForward<decltype(evalFunc), AdjointData*&>;

      AdjointInterfaceImpl<Real, Index, AdjointData> interface(adjointData);

      evaluateExtFuncForward(start, end, forwardFunc, jacobiVector, &interface, evalFunc, adjointData);
    }

  public:

    /**
     * @brief Register a variable as an active variable.
     *
     * The index of the variable is set to a non zero number.
     * @param[in,out] value The value which will be marked as an active variable.
     */
    CODI_INLINE void registerInput(ActiveReal<JacobiIndexTape<TapeTypes> >& value) {
      indexHandler.assignUnusedIndex(value.getGradientData());
    }

    /**
     * @brief Modify the output of an external function such that the tape sees it as an active variable.
     *
     * @param[in,out] value  The output value of the external function.
     * @return Zero
     */
    CODI_INLINE Real registerExtFunctionOutput(ActiveReal<JacobiIndexTape<TapeTypes> >& value) {
      registerInput(value);

      return Real();
    }

    /**
     * @brief Not needed in this implementation.
     *
     * @param[in] value Not used.
     */
    CODI_INLINE void registerOutput(ActiveReal<JacobiIndexTape<TapeTypes> >& value) {
      if(!IndexHandler::AssignNeedsStatement) {
        value = PassiveReal(1.0) * value;
      }
    }

    /**
     * @brief Gather the general performance values of the tape.
     *
     * Computes values like the total amount of statements stored or the currently consumed memory.
     *
     * @return The values for the tape.
     */
    TapeValues getTapeValues() const {
      std::string name = "CoDi Tape Statistics (" + std::string(TapeTypes::tapeName) + ")";
      TapeValues values(name);

      addTapeBaseValues(values);
      addStmtValues(values);
      addJacobiValues(values);
      addExtFuncValues(values);

      return values;
    }
  };

  /**
   * @brief The instantiation of the index manager for the jacobi index tapes.
   */
  template <typename TapeTypes>
  typename TapeTypes::IndexHandler JacobiIndexTape<TapeTypes>::indexHandler(0);

}
