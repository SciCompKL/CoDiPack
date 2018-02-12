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

#include <cstddef>
#include <tuple>

#include "../activeReal.hpp"
#include "../expressionHandle.hpp"
#include "chunkVector.hpp"
#include "indices/reuseIndexHandler.hpp"
#include "handles/functionHandleFactory.hpp"
#include "handles/staticObjectHandleFactory.hpp"
#include "handles/staticFunctionHandleFactory.hpp"
#include "primalTapeExpressions.hpp"
#include "reverseTapeInterface.hpp"
#include "singleChunkVector.hpp"
#include "../tools/tapeValues.hpp"

namespace codi {

  /**
   * @brief Vector definition for the ChunkPrimalValueIndexTape.
   *
   * The structure defines all vectors as chunk vectors.
   *
   * See PrimalValueIndexTape for details.
   *
   * @tparam          Real  The type for the primal values.
   * @tparam  IndexHandler  The index handler for the managing of the indices. It has to be a index handler that assumes index reuse.
   * @tparam GradientValue  The type for the adjoint values. (Default: Same as the primal value.)
   * @tparam HandleFactory  The factory for the reverse interpretation of the expressions. Needs to implement the HandleFactoryInterface class.
   * @tparam    DataVector  The data manager for the chunks. Needs to implement a ChunkVector interface.
   */
  template <typename RTT, template<typename> class HandleFactoryType, template<typename, typename> class DataVector>
  struct IndexPrimalValueTapeTypes : public RTT {

    CODI_INLINE_REVERSE_TAPE_TYPES(RTT)

    /** @brief The factory for the expression handles. */
    typedef HandleFactoryType<RTT> HandleFactory;

    /** @brief The data type for the created handles. */
    typedef typename HandleFactory::Handle Handle;

    /** @brief The tape type structure, that defines the basic types. */
    typedef RTT BaseTypes;

    /** @brief The data for each statement. */
    typedef Chunk4<Index, Real, Handle, StatementInt> StatementChunk;
    /** @brief The chunk vector for the statement data. */
    typedef DataVector<StatementChunk, EmptyChunkVector> StatementVector;

    /** @brief The data for the indices of each statement */
    typedef Chunk1< Index> IndexChunk;
    /** @brief The chunk vector for the index data. */
    typedef DataVector<IndexChunk, StatementVector> IndexVector;

    /** @brief The data for the constant values of each statement */
    typedef Chunk1< PassiveReal> ConstantValueChunk;
    /** @brief The chunk vector for the constant data. */
    typedef DataVector<ConstantValueChunk, IndexVector> ConstantValueVector;

    /** @brief The data for the external functions. */
    typedef Chunk2<ExternalFunction,typename ConstantValueVector::Position> ExternalFunctionChunk;
    /** @brief The chunk vector for the external  function data. */
    typedef DataVector<ExternalFunctionChunk, ConstantValueVector> ExternalFunctionVector;

    /** @brief The position for all the different data vectors. */
    typedef typename ExternalFunctionVector::Position Position;

    /** @brief The name of the tape as a string. */
    constexpr static const char* tapeName = "PrimalValueIndexTape";

  };

  /**
   * @brief A reverse AD tape that stores primal values for the reverse evaluation.
   *
   * The PrimalValueIndexTape implements a fully featured ReverseTapeInterface. Depending on
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
  class PrimalValueIndexTape : public ReverseTapeInterface<typename TapeTypes::Real, typename TapeTypes::Index, typename TapeTypes::GradientValue, PrimalValueIndexTape<TapeTypes>, typename TapeTypes::Position > {
  public:

    CODI_INLINE_REVERSE_TAPE_TYPES(TapeTypes::BaseTypes)

    /** @brief The factory for the expression handles. */
    typedef typename TapeTypes::HandleFactory HandleFactory;

    /** @brief The data type for the created handles. */
    typedef typename HandleFactory::Handle Handle;

    /** @brief The tape type structure, that defines the basic types. */
    typedef typename TapeTypes::BaseTypes BaseTypes;

    /** @brief The gradient data is just the index type. */
    typedef Index GradientData;

    /** @brief The termination of the vector sequence. */
    EmptyChunkVector emptyVector;

    /** @brief The index handler for the active real's. */
    static typename TapeTypes::IndexHandler indexHandler;

    /** @brief Disables code path in CoDiPack that are optimized for Jacobi taping */
    static const bool AllowJacobiOptimization = false;

    /** @brief This tape requires no special primal value handling since the primal value vector is not overwritten. */
    static const bool RequiresPrimalReset = true;

    #define TAPE_NAME PrimalValueIndexTape

    #define POSITION_TYPE typename TapeTypes::Position
    #define INDEX_HANDLER_NAME indexHandler
    #define RESET_FUNCTION_NAME resetPrimalValues
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

    #define ROOT_VECTOR extFuncVector
    #include "modules/ioModule.tpp"

    // TAPE_NAME is undefined at the end of the file

    /** @brief The temporary vector for the reverse evaluation. */
    Real* primalsCopy;

    /** @brief The size of the copied primal vector */
    Index primalsCopySize;

  public:
    /**
     * @brief Creates a tape with the size of zero for the data, statements and external functions.
     */
    PrimalValueIndexTape() :
      emptyVector(),
      /* defined in tapeBaseModule */adjoints(NULL),
      /* defined in tapeBaseModule */adjointsSize(0),
      /* defined in tapeBaseModule */active(false),
      /* defined in the primalValueModule */stmtVector(DefaultChunkSize, &emptyVector),
      /* defined in the primalValueModule */indexVector(DefaultChunkSize, &stmtVector),
      /* defined in the primalValueModule */constantValueVector(DefaultChunkSize, &indexVector),
      /* defined in the primalValueModule */primals(NULL),
      /* defined in the primalValueModule */primalsSize(0),
      /* defined in the primalValueModule */primalsIncr(DefaultSmallChunkSize),
      /* defined in externalFunctionsModule */extFuncVector(1000, &constantValueVector),
      primalsCopy(NULL),
      primalsCopySize(0) {}

    /** @brief Tear down the tape. Delete all values from the modules */
    ~PrimalValueIndexTape() {
      cleanTapeBase();

      if(NULL != primalsCopy) {
        free(primalsCopy);
        primalsCopy = NULL;
      }

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
    void swap(PrimalValueIndexTape& other) {
      swapTapeBaseModule(other);
      swapPrimalValueModule(other);

      // the index handler is not swapped because the indices of the program state need to stay valid

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
      CODI_UNUSED(start);
      CODI_UNUSED(end);
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
    CODI_INLINE void pushStmtData(Index& lhsIndex, const Real& rhsValue, const Handle& handle, const StatementInt& passiveVariableNumber) {
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
    CODI_INLINE void store(Real& lhsValue, Index& lhsIndex, const ActiveReal<PrimalValueIndexTape<TapeTypes> >& rhs) {
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
     * @brief Evaluate the stack from the start to to the end position for the primal evaluation.
     *
     * It has to hold start <= end.
     *
     * @param[in,out]       stmtPos  The current position in the statement data. It will decremented in the method.
     * @param[in]        endStmtPos  The ending position for statement data.
     * @param[in]        lhsIndices  The indices from the lhs of each statement.
     * @param[in]     storedPrimals  The overwritten primal from the primal vector.
     * @param[in]        statements  The vector with the handles for each statement.
     * @param[in] passiveActiveReal  The number passive values for each statement.
     * @param[in,out]      indexPos  The current position for the index data. It will decremented in the method.
     * @param[in]           indices  The indices for the arguments of the rhs.
     * @param[in,out]   constantPos  The current position in the constant data vector. It will decremented in the method.
     * @param[in]         constants  The constant values in the rhs expressions.
     * @param[in,out]  primalVector  The vector of the primal variables.
     *
     * @tparam AdjointData The data for the adjoint vector it needs to support add, multiply and comparison operations.
     */
    CODI_INLINE void evaluateStackPrimal(Real* primalVector,
                                         size_t& constantPos, const size_t& endConstantPos, PassiveReal* &constants,
                                         size_t& indexPos, const size_t& endIndexPos, Index* &indices,
                                         size_t& stmtPos, const size_t& endStmtPos, Index* lhsIndices, Real* storedPrimals,
                                         Handle* &statements, StatementInt* &passiveActiveReal) {
      CODI_UNUSED(endConstantPos);
      CODI_UNUSED(endIndexPos);

      while(stmtPos < endStmtPos) {
        const Index& lhsIndex = lhsIndices[stmtPos];

        storedPrimals[stmtPos] = primalVector[lhsIndex];
        primalVector[lhsIndex] = HandleFactory::template callPrimalHandle<PrimalValueIndexTape<TapeTypes> >(statements[stmtPos], passiveActiveReal[stmtPos], indexPos, indices, constantPos, constants, primalVector);
        stmtPos += 1;
      }
    }

    /**
     * @brief Evaluate the stack from the start to to the end position.
     *
     * It has to hold start >= end.
     *
     * @param[in,out]       stmtPos  The current position in the statement data. It will decremented in the method.
     * @param[in]        endStmtPos  The ending position for statement data.
     * @param[in]        lhsIndices  The indices from the lhs of each statement.
     * @param[in]     storedPrimals  The overwritten primal from the primal vector.
     * @param[in]        statements  The vector with the handles for each statement.
     * @param[in] passiveActiveReal  The number passive values for each statement.
     * @param[in,out]      indexPos  The current position for the index data. It will decremented in the method.
     * @param[in]           indices  The indices for the arguments of the rhs.
     * @param[in,out]   constantPos  The current position in the constant data vector. It will decremented in the method.
     * @param[in]         constants  The constant values in the rhs expressions.
     * @param[in,out] adjointData  The vector of the adjoint varaibles.
     *
     * @tparam AdjointData The data for the adjoint vector it needs to support add, multiply and comparison operations.
     */
    template<typename AdjointData>
    CODI_INLINE void evaluateStack(AdjointData* adjointData, Real* primalVector,
                                   size_t& constantPos, const size_t& endConstantPos, PassiveReal* &constants,
                                   size_t& indexPos, const size_t& endIndexPos, Index* &indices,
                                   size_t& stmtPos, const size_t& endStmtPos, Index* lhsIndices, Real* storedPrimals,
                                      Handle* &statements, StatementInt* &passiveActiveReal) {
      CODI_UNUSED(endConstantPos);
      CODI_UNUSED(endIndexPos);

      while(stmtPos > endStmtPos) {
        --stmtPos;
        const Index& lhsIndex = lhsIndices[stmtPos];

        primalVector[lhsIndex] = storedPrimals[stmtPos];
#if CODI_EnableVariableAdjointInterfaceInPrimalTapes
          adjointData->setLhsAdjoint(lhsIndex);
          adjointData->resetAdjointVec(lhsIndex);

          HandleFactory::template callHandle<PrimalValueIndexTape<TapeTypes> >(statements[stmtPos], 1.0, passiveActiveReal[stmtPos], indexPos, indices, constantPos, constants, primalVector, adjointData);
#else
          const GradientValue adj = adjointData[lhsIndex];
          adjointData[lhsIndex] = GradientValue();

          HandleFactory::template callHandle<PrimalValueIndexTape<TapeTypes> >(statements[stmtPos], adj, passiveActiveReal[stmtPos], indexPos, indices, constantPos, constants, primalVector, adjointData);
#endif
      }
    }

    /**
     * @brief Evaluate a part of the statement vector for the primal.
     *
     * It has to hold start <= end.
     *
     * The function calls the primal  evaluation method for the stack.
     *
     * @param[in] start  The starting point for the statement vector.
     * @param[in]   end  The ending point for the statement vector.
     * @param[in]  args  The arguments from the other vectors.
     *
     * @tparam Args  The types of the other arguments.
     */
    template<typename ... Args>
    CODI_INLINE void evalStmtPrimal(const StmtPosition& start, const StmtPosition& end, Args&&... args) {
      Index* data1;
      Real* data2;
      Handle* data3;
      StatementInt* data4;

      size_t dataPos = start.data;
      for(size_t curChunk = start.chunk; curChunk < end.chunk; ++curChunk) {
        stmtVector.getDataAtPosition(curChunk, 0, data1, data2, data3, data4);

        evaluateStackPrimal(dataPos, stmtVector.getChunkUsedData(curChunk), data1, data2, data3, data4, std::forward<Args>(args)...);

        codiAssert(dataPos == stmtVector.getChunkUsedData(curChunk)); // After a full chunk is evaluated the data position needs to be at the end

        dataPos = 0;
      }

      // Iterate over the reminder also covers the case if the start chunk and end chunk are the same
      stmtVector.getDataAtPosition(end.chunk, 0, data1, data2, data3, data4);
      evaluateStackPrimal(dataPos, end.data, data1, data2, data3, data4, std::forward<Args>(args)...);
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
      Index* data1;
      Real* data2;
      Handle* data3;
      StatementInt* data4;

      size_t dataPos = start.data;
      for(size_t curChunk = start.chunk; curChunk > end.chunk; --curChunk) {
        stmtVector.getDataAtPosition(curChunk, 0, data1, data2, data3, data4);

        evaluateStack(dataPos, 0, data1, data2, data3, data4, std::forward<Args>(args)...);

        codiAssert(dataPos == 0); // after a full chunk is evaluated, the data position needs to be zero

        dataPos = stmtVector.getChunkUsedData(curChunk - 1);
      }

      // Iterate over the reminder also covers the case if the start chunk and end chunk are the same
      stmtVector.getDataAtPosition(end.chunk, 0, data1, data2, data3, data4);
      evaluateStack(dataPos, end.data, data1, data2, data3, data4, std::forward<Args>(args)...);
    }

    /**
     * @brief Evaluate a part of the statement vector for the primal.
     *
     * It has to hold start <= end.
     *
     * The function calls the primal evaluation method for the constant value vector.
     *
     * @param[in] start  The starting point for the statement vector.
     * @param[in]   end  The ending point for the statement vector.
     * @param[in]  args  The arguments from the other vectors.
     *
     * @tparam Args  The types of the other arguments.
     */
    template<typename ... Args>
    CODI_INLINE void evalExtFuncPrimalCallback(const ConstantValuePosition& start, const ConstantValuePosition& end, Args&&... args) {
      evaluateConstantValuesPrimal(start, end, std::forward<Args>(args)...);
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

    /**
     * @brief Allocates a copy of the primal vector that is used in the evaluation.
     *
     * It has to hold start >= end.
     *
     * The function calls the evaluation method for the external function vector.
     *
     * @param[in]            start  The starting point for the statement vector.
     * @param[in]              end  The ending point for the statement vector.
     * @param[in,out]  adjointData  The adjoint vector for the evaluation.
     *
     * @tparam AdjointData  The data type for the adjoint vector.
     */
    template<typename AdjointData>
    CODI_INLINE void evaluateInt(const Position& start, const Position& end, AdjointData* adjointData) {
      if(primalsCopySize < primalsSize) {
        primalsCopy = (Real*)realloc(primalsCopy, sizeof(Real) * primalsSize);
        primalsCopySize = primalsSize;
      }
      memcpy(primalsCopy, primals, sizeof(Real) * primalsSize);

      AdjointInterfacePrimalImpl<Real, AdjointData> interface(adjointData, primalsCopy);


#if CODI_EnableVariableAdjointInterfaceInPrimalTapes
      typedef AdjointInterfaceImpl<Real, AdjointData> AdjVecType;
      AdjVecType* adjVec = &interface;
#else
      static_assert(std::is_same<AdjointData, GradientValue>::value,
        "Please enable 'CODI_EnableVariableAdjointInterfaceInPrimalTapes' in order"
        " to use custom adjoint vectors in the primal value tapes.");

      typedef AdjointData AdjVecType;
      AdjVecType* adjVec = adjointData;
#endif

      auto evalFunc = [this] (AdjVecType* adjointData, Real* primalVector,
                              size_t& constantPos, const size_t& endConstantPos, PassiveReal* &constants,
                              size_t& indexPos, const size_t& endIndexPos, Index* &indices,
                              size_t& stmtPos, const size_t& endStmtPos, Index* lhsIndices, Real* storedPrimals,
                                Handle* &statements, StatementInt* &passiveActiveReal) {
        evaluateStack<AdjVecType>(adjointData, primalVector, constantPos, endConstantPos, constants,
                                     indexPos, endIndexPos, indices, stmtPos, endStmtPos, lhsIndices, storedPrimals,
                                     statements, passiveActiveReal);
      };
      auto reverseFunc = &ConstantValueVector::template evaluateReverse<decltype(evalFunc), AdjVecType*&, Real*&>;
      evaluateExtFunc(start, end, reverseFunc, constantValueVector, &interface, evalFunc, adjVec, primalsCopy);
    }

    /**
     * Function object that is used to revert the primal values in the primal value vector.
     */
    struct PrimalValueReseter {
      PrimalValueIndexTape& tape;

      /**
       * @brief Create the function object.
       *
       * @param[in,out]     tape  The reference to the actual tape.
       */
      PrimalValueReseter(PrimalValueIndexTape& tape) :
        tape(tape){}

      /**
       * @brief The operator reverts the value in the primal value vector to the old value stored during the recording.
       *
       * @param[in]    index  The lhs index for the statement.
       * @param[in]    value  The old value that has been overwritten in the primal value vector.
       * @param[in]   handle  The handle for the evaluation of the statement.
       * @param[in] stmtSize  The number of non zero entries in the statement.
       */
      void operator () (Index* index, Real* value, Handle* handle, StatementInt* stmtSize) {
        CODI_UNUSED(handle);
        CODI_UNUSED(stmtSize);

        tape.primals[*index] = *value;
      }
    };


    /**
     * @brief Resets the primal values up to the specified position.
     *
     * The stack is reversed such that the primal value changes are reversed until the specified position is reached.
     * The primal value vector of the tape has the same status when the position was recorded.
     *
     * @param[in] pos  The position for the tape reset.
     */
    CODI_INLINE void resetPrimalValues(const Position& pos) {

      // Do not perform a global reset on the primal value vector if the tape is cleared
      if(getZeroPosition() != pos) {

        PrimalValueReseter reseter(*this);

        StmtPosition stmtEnd = stmtVector.getPosition();

        stmtVector.forEach(stmtEnd, pos.inner.inner.inner, reseter);
      }

      // call the function from the external function module
      resetExtFunc(pos);
    }

  public:


    /**
     * @brief Specialized evaluate function for the preaccumulation.
     *
     * The function does a normal reverse evaluation and afterwards the primal value vector es reset to the
     * current state.
     *
     * It has to hold start >= end
     *
     * @param[in] start  The starting position for the evaluation.
     * @param[in]   end  The stopping position for the evaluation.
     */
    CODI_INLINE void evaluatePreacc(const Position& start, const Position& end) {

      resizeAdjointsToIndexSize();

      std::swap(primals, primalsCopy);

      AdjointInterfacePrimalImpl<Real, GradientValue> interface(adjoints, primalsCopy);

#if CODI_EnableVariableAdjointInterfaceInPrimalTapes
        typedef AdjointInterfaceImpl<Real, GradientValue> AdjVecType;
        AdjVecType* adjVec = &interface;
#else
        typedef GradientValue AdjVecType;
        AdjVecType* adjVec = adjoints;
#endif
      auto reverseFunc = [this] (AdjVecType* adjointData, Real* primalVector,
                                 size_t& constantPos, const size_t& endConstantPos, PassiveReal* &constants,
                                 size_t& indexPos, const size_t& endIndexPos, Index* &indices,
                                 size_t& stmtPos, const size_t& endStmtPos, Index* lhsIndices, Real* storedPrimals,
                                   Handle* &statements, StatementInt* &passiveActiveReal) {
        evaluateStack<AdjVecType>(adjointData, primalVector, constantPos, endConstantPos, constants,
                                     indexPos, endIndexPos, indices, stmtPos, endStmtPos, lhsIndices, storedPrimals,
                                     statements, passiveActiveReal);
      };
      auto reverseIter = &ConstantValueVector::template evaluateReverse<decltype(reverseFunc), AdjVecType*&, Real*&>;
      evaluateExtFunc(start, end, reverseIter, constantValueVector, &interface, reverseFunc, adjVec, primalsCopy);

      auto primalFunc = [this] (Real* primalVector,
                                size_t& constantPos, const size_t& endConstantPos, PassiveReal* &constants,
                                size_t& indexPos, const size_t& endIndexPos, Index* &indices,
                                size_t& stmtPos, const size_t& endStmtPos, Index* lhsIndices, Real* storedPrimals,
                                  Handle* &statements, StatementInt* &passiveActiveReal) {
        evaluateStackPrimal(primalVector, constantPos, endConstantPos, constants,
                            indexPos, endIndexPos, indices, stmtPos, endStmtPos, lhsIndices, storedPrimals,
                            statements, passiveActiveReal);
      };
      auto primalIter = &ConstantValueVector::template evaluateForward<decltype(primalFunc), Real*&>;
      evaluateExtFuncPrimal(end, start, primalIter, constantValueVector, primalFunc, primalsCopy);

      std::swap(primals, primalsCopy);
    }

    /**
     * @brief Register a variable as an active variable.
     *
     * The index of the variable is set to a non zero number.
     * @param[in,out] value The value which will be marked as an active variable.
     */
    CODI_INLINE void registerInput(ActiveReal<PrimalValueIndexTape<TapeTypes> >& value) {
      if(isActive()) {
        //TODO: Test if register of already active value will yield errors.
        indexHandler.assignIndex(value.getGradientData());

        checkPrimalsSize();
        primals[value.getGradientData()] = value.getValue();
      }
    }


    /**
     * @brief Modify the output of an external function such that the tape sees it as an active variable.
     *
     * @param[in,out] value  The output value of the external function.
     * @return The old internally stored value for the newly generated index.
     */
    CODI_INLINE Real registerExtFunctionOutput(ActiveReal<PrimalValueIndexTape<TapeTypes> >& value) {
      Real oldValue = value.getValue();
      if(isActive()) {
        //TODO: Test if register of already active value will yield errors.
        indexHandler.assignIndex(value.getGradientData());

        checkPrimalsSize();
        oldValue = primals[value.getGradientData()];
        primals[value.getGradientData()] = value.getValue();
      }

      return oldValue;
    }

    /**
     * @brief Helper function for external functions.
     *
     * The function sets a primal value in current primal value vector. The value that needs to be provided here is the
     * one that is returned from registerExtFunctionOutput.
     *
     * @param[in]  index  The index on which the primal value is changed.
     * @param[in] primal  The new primal value for the index.
     */
    CODI_INLINE void setExternalValueChange(const GradientData& index, const Real& primal) {

      // TODO: only allowed during reverse evaluation. Add check.

      primalsCopy[index] = primal;
    }

    /**
     * @brief It is ensured that each output variable has a unique index.
     *
     * @param[in] value  The value will have an unique index that is used by no other variable.
     */
    CODI_INLINE void registerOutput(ActiveReal<PrimalValueIndexTape<TapeTypes> >& value) {
      if(isActive() && value.getGradientData() != 0) {
        if(!IndexHandler::AssignNeedsStatement) {
          Index rhsIndex = value.getGradientData();

          indexHandler.assignIndex(value.getGradientData());
          pushCopyHandle(value.getValue(), value.getGradientData(), rhsIndex);
        }
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
      addPrimalValueValues(values);
      addExtFuncValues(values);

      return values;
    }
  };

  /**
   * @brief The instantiation of the index manager for the primal index tapes.
   */
  template <typename TapeTypes>
  typename TapeTypes::IndexHandler PrimalValueIndexTape<TapeTypes>::indexHandler(MaxStatementIntSize - 1);

  #include "modules/primalValueStaticModule.tpp"
  #undef TAPE_NAME

}
