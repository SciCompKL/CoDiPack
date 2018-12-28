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
#include "indices/linearIndexHandler.hpp"
#include "handles/functionHandleFactory.hpp"
#include "primalTapeExpressions.hpp"
#include "reverseTapeInterface.hpp"
#include "singleChunkVector.hpp"
#include "../tapeTypes.hpp"
#include "../tools/tapeValues.hpp"

namespace codi {

  /**
   * @brief Vector definition for the ChunkPrimalValueTape.
   *
   * The structure defines all vectors as chunk vectors.
   *
   * See PrimalValueTape for details.
   *
   * @tparam           RTT  The basic type definitions for the tape. Need to define everything from ReverseTapeTypes.
   * @tparam HandleFactory  The factory for the reverse interpretation of the expressions. Needs to implement the HandleFactoryInterface class.
   * @tparam    DataVector  The data manager for the chunks. Needs to implement a ChunkVector interface.
   */
  template <typename RTT, template<typename> class HandleFactoryType, template<typename, typename> class DataVector>
  struct PrimalValueTapeTypes {

    CODI_INLINE_REVERSE_TAPE_TYPES(RTT)

    /** @brief The factory for the expression handles. */
    typedef HandleFactoryType<RTT> HandleFactory;

    /** @brief The data type for the created handles. */
    typedef typename HandleFactory::Handle Handle;

    /** @brief The tape type structure, that defines the basic types. */
    typedef RTT BaseTypes;

    /** @brief The data for each statement. */
    typedef Chunk2<Handle, StatementInt> StatementChunk;
      /** @brief The chunk vector for the statement data. */
    typedef DataVector<StatementChunk, IndexHandler> StatementVector;

    /** @brief The data for the indices of each statement */
    typedef Chunk1< Index> IndexChunk;
    /** @brief The chunk vector for the index data. */
    typedef DataVector<IndexChunk, StatementVector> IndexVector;

    /** @brief The data for the constant values of each statement */
    typedef Chunk1< Real> PassiveValueChunk;
    /** @brief The chunk vector for the constant data. */
    typedef DataVector<PassiveValueChunk, IndexVector> PassiveValueVector;

    /** @brief The data for the constant values of each statement */
    typedef Chunk1< PassiveReal> ConstantValueChunk;
    /** @brief The chunk vector for the constant data. */
    typedef DataVector<ConstantValueChunk, PassiveValueVector> ConstantValueVector;

    /** @brief The data for the external functions. */
    typedef Chunk2<ExternalFunction,typename ConstantValueVector::Position> ExternalFunctionChunk;
    /** @brief The chunk vector for the external  function data. */
    typedef DataVector<ExternalFunctionChunk, ConstantValueVector> ExternalFunctionVector;

    /** @brief The position for all the different data vectors. */
    typedef typename ExternalFunctionVector::Position Position;

    /** @brief The name of the tape as a string. */
    constexpr static const char* tapeName = "PrimalValueTape";

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
  class PrimalValueTape : public ReverseTapeInterface<typename TapeTypes::Real, typename TapeTypes::Index, typename TapeTypes::GradientValue, PrimalValueTape<TapeTypes>, typename TapeTypes::Position >  {
  public:

    CODI_INLINE_REVERSE_TAPE_TYPES(TapeTypes::BaseTypes)

    /** @brief The tape type structure, that defines the basic types. */
    typedef typename TapeTypes::BaseTypes BaseTypes;

    /** @brief The generator for all handles to expression templates. */
    typedef typename TapeTypes::HandleFactory HandleFactory;

    /** @brief The actual type of the stored handles. */
    typedef typename HandleFactory::Handle Handle;

    /** @brief The gradient data is just the index type. */
    typedef Index GradientData;

    /** @brief The index handler for the active real's. */
    IndexHandler indexHandler;

    /** @brief Disables code path in CoDiPack that are optimized for Jacobi taping */
    static const bool AllowJacobiOptimization = false;

    /** @brief This tape requires no special primal value handling since the primal value vector is not overwritten. */
    static const bool RequiresPrimalReset = false;

    #define TAPE_NAME PrimalValueTape

    #define POSITION_TYPE typename TapeTypes::Position
    #define INDEX_HANDLER_NAME indexHandler
    #define RESET_FUNCTION_NAME resetExtFunc
    #define EVALUATE_FUNCTION_NAME evaluateInt
    #define EVALUATE_FORWARD_FUNCTION_NAME evaluateForwardInt
    #define EVALUATE_PRIMAL_FUNCTION_NAME evaluatePrimalInt
    #include "modules/tapeBaseModule.tpp"

    #define CHILD_VECTOR_TYPE EmptyChunkVector
    #define STMT_VECTOR_TYPE typename TapeTypes::StatementVector
    #define INDEX_VECTOR_TYPE typename TapeTypes::IndexVector
    #define PASSIVE_VECTOR_TYPE typename TapeTypes::PassiveValueVector
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
      /* defined in the primalValueModule */passiveValueVector(DefaultChunkSize, &indexVector),
      /* defined in the primalValueModule */constantValueVector(DefaultChunkSize, &passiveValueVector),
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

      Index startPos = min(end.inner.inner.inner.inner.inner, adjointsSize - 1);
      Index endPos = min(start.inner.inner.inner.inner.inner, adjointsSize - 1);

      for(Index i = startPos + 1; i <= endPos; ++i) {
        adjoints[i] = GradientValue();
      }
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

      AdjVecInterface<AdjointData> interface(adjointData, primals);
      AdjVecType<AdjointData>* adjVec = wrapAdjointVector(interface, adjointData);

      auto evalFunc = [this] (const size_t& startAdjPos, const size_t& endAdjPos,
                              AdjVecType<AdjointData>* adjointData,
                              size_t& constantPos, const size_t& endConstantPos, PassiveReal* &constants,
                              size_t& passivePos, const size_t& endPassivePos, Real* &passives,
                              size_t& indexPos, const size_t& endIndexPos, Index* &indices,
                              size_t& stmtPos, const size_t& endStmtPos,
                                Handle* &statements, StatementInt* &passiveActiveReal) {
        evaluateStackReverse<AdjVecType<AdjointData>>(startAdjPos, endAdjPos, adjointData, constantPos, endConstantPos, constants, passivePos, endPassivePos, passives,
                                     indexPos, endIndexPos, indices, stmtPos, endStmtPos, statements, passiveActiveReal);
      };
      auto reverseFunc = &ConstantValueVector::template evaluateReverse<decltype(evalFunc), AdjVecType<AdjointData>*&>;
      evaluateExtFunc(start, end, reverseFunc, constantValueVector, &interface, evalFunc, adjVec);
    }

    /**
     * @brief Allocates a copy of the primal vector that is used in the evaluation.
     *
     * It has to hold start <= end.
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
    CODI_INLINE void evaluateForwardInt(const Position& start, const Position& end, AdjointData* adjointData) {

      AdjVecInterface<AdjointData> interface(adjointData, primals);
      AdjVecType<AdjointData>* adjVec = wrapAdjointVector(interface, adjointData);

      auto evalFunc = [this] (const size_t& startAdjPos, const size_t& endAdjPos,
                              AdjVecType<AdjointData>* adjointData,
                              size_t& constantPos, const size_t& endConstantPos, PassiveReal* &constants,
                              size_t& passivePos, const size_t& endPassivePos, Real* &passives,
                              size_t& indexPos, const size_t& endIndexPos, Index* &indices,
                              size_t& stmtPos, const size_t& endStmtPos,
                                Handle* &statements, StatementInt* &passiveActiveReal) {
        evaluateStackForward<AdjVecType<AdjointData>>(startAdjPos, endAdjPos, adjointData, constantPos, endConstantPos, constants, passivePos, endPassivePos, passives,
                                     indexPos, endIndexPos, indices, stmtPos, endStmtPos, statements, passiveActiveReal);
      };
      auto forwardFunc = &ConstantValueVector::template evaluateForward<decltype(evalFunc), AdjVecType<AdjointData>*&>;
      evaluateExtFuncForward(start, end, forwardFunc, constantValueVector, &interface, evalFunc, adjVec);
    }

    /**
     * @brief Evaluate the tape from start to end.
     *
     * The function performs the primal evaluation of the recorded tape from
     * the start position to the end position.
     *
     * The primal evaluation will update the internal primal value vector.
     *
     * It has to hold start <= end.
     *
     * @param[in] start The starting position for the forward evaluation.
     * @param[in]   end The ending position for the forward evaluation.
     */
    CODI_INLINE void evaluatePrimalInt(const Position& start, const Position& end) {

      AdjVecInterface<GradientValue> interface(adjoints, primals);

      auto evalFunc = [this] (const size_t& startAdjPos, const size_t& endAdjPos,
                              size_t& constantPos, const size_t& endConstantPos, PassiveReal* &constants,
                              size_t& passivePos, const size_t& endPassivePos, Real* &passives,
                              size_t& indexPos, const size_t& endIndexPos, Index* &indices,
                              size_t& stmtPos, const size_t& endStmtPos, Handle* &statements, StatementInt* &passiveActiveReal) {
        evaluateStackPrimal(startAdjPos, endAdjPos, constantPos, endConstantPos, constants, passivePos, endPassivePos, passives,
                            indexPos, endIndexPos, indices, stmtPos, endStmtPos, statements, passiveActiveReal);
      };
      auto primalFunc = &ConstantValueVector::template evaluateForward<decltype(evalFunc)>;
      evaluateExtFuncPrimal(start, end, primalFunc, constantValueVector, &interface, evalFunc);
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
    CODI_INLINE void store(Real& lhsValue, Index& lhsIndex, const ActiveReal<PrimalValueTape<TapeTypes> >& rhs) {
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
     * @param[in,out]   adjointData  The vector of the adjoint variables.
     * @param[in,out]   constantPos  The current position in the constant data vector. It will decremented in the method.
     * @param[in]    endConstantPos  The ending position in the constant data vector.
     * @param[in]         constants  The constant values in the rhs expressions.
     * @param[in,out]      indexPos  The current position for the index data. It will decremented in the method.
     * @param[in]       endIndexPos  The ending position in the index data.
     * @param[in]           indices  The indices for the arguments of the rhs.
     * @param[in,out]       stmtPos  The current position in the statement data. It will decremented in the method.
     * @param[in]        statements  The vector with the handles for each statement.
     * @param[in] passiveActiveReal  The number passive values for each statement.
     *
     * @tparam AdjointData The data for the adjoint vector it needs to support add, multiply and comparison operations.
     */
    template<typename AdjointData>
    CODI_INLINE void evaluateStackReverse(const size_t& startAdjPos, const size_t& endAdjPos, AdjointData* adjointData,
                                          size_t& constantPos, const size_t& endConstPos, PassiveReal* &constants,
                                          size_t& passivePos, const size_t& endPassivePos, Real* &passives,
                                          size_t& indexPos, const size_t& endIndexPos, Index* &indices,
                                          size_t& stmtPos, const size_t& endStmtPos, Handle* &statements,
                                          StatementInt* &passiveActiveReal) {
      CODI_UNUSED(endConstPos);
      CODI_UNUSED(endPassivePos);
      CODI_UNUSED(endIndexPos);
      CODI_UNUSED(endStmtPos);

      size_t adjPos = startAdjPos;

      while(adjPos > endAdjPos) {
        --stmtPos;
#if CODI_EnableVariableAdjointInterfaceInPrimalTapes
          adjointData->setLhsAdjoint(adjPos);
          if(ZeroAdjointReverse && StatementIntInputTag != passiveActiveReal[stmtPos]) {
            adjointData->resetAdjointVec(adjPos);
          }
#else
          const AdjointData adj = adjointData[adjPos];
          if(ZeroAdjointReverse && StatementIntInputTag != passiveActiveReal[stmtPos]) {
            adjointData[adjPos] = GradientValue();
          }
#endif
        --adjPos;

        if(StatementIntInputTag != passiveActiveReal[stmtPos]) {
#if CODI_EnableVariableAdjointInterfaceInPrimalTapes
          HandleFactory::template callHandle<PrimalValueTape<TapeTypes> >(statements[stmtPos], 1.0, passiveActiveReal[stmtPos], indexPos, indices, passivePos, passives, constantPos, constants, primals, adjointData);
#else
          HandleFactory::template callHandle<PrimalValueTape<TapeTypes> >(statements[stmtPos], adj, passiveActiveReal[stmtPos], indexPos, indices, passivePos, passives, constantPos, constants, primals, adjointData);
#endif
        }
      }
    }

    /**
     * @brief Evaluate the stack from the start to to the end position.
     *
     * It has to hold start <= end.
     *
     * @param[in]       startAdjPos  The starting position for the adjoint evaluation.
     * @param[in]         endAdjPos  The ending position for the adjoint evaluation.
     * @param[in,out]   adjointData  The vector of the adjoint variables.
     * @param[in,out]   constantPos  The current position in the constant data vector. It will decremented in the method.
     * @param[in]    endConstantPos  The ending position in the constant data vector.
     * @param[in]         constants  The constant values in the rhs expressions.
     * @param[in,out]      indexPos  The current position for the index data. It will decremented in the method.
     * @param[in]       endIndexPos  The ending position in the index data.
     * @param[in]           indices  The indices for the arguments of the rhs.
     * @param[in,out]       stmtPos  The current position in the statement data. It will decremented in the method.
     * @param[in]        statements  The vector with the handles for each statement.
     * @param[in] passiveActiveReal  The number passive values for each statement.
     *
     * @tparam AdjointData The data for the adjoint vector it needs to support add, multiply and comparison operations.
     */
    template<typename AdjointData>
    CODI_INLINE void evaluateStackForward(const size_t& startAdjPos, const size_t& endAdjPos, AdjointData* adjointData,
                                          size_t& constantPos, const size_t& endConstPos, PassiveReal* &constants,
                                         size_t& passivePos, const size_t& endPassivePos, Real* &passives,
                                          size_t& indexPos, const size_t& endIndexPos, Index* &indices,
                                          size_t& stmtPos, const size_t& endStmtPos, Handle* &statements,
                                          StatementInt* &passiveActiveReal) {
      CODI_UNUSED(endConstPos);
      CODI_UNUSED(endPassivePos);
      CODI_UNUSED(endIndexPos);
      CODI_UNUSED(endStmtPos);

      size_t adjPos = startAdjPos;

      while(adjPos < endAdjPos) {
        adjPos += 1;

        GradientValue lhsAdj = GradientValue();

        if(StatementIntInputTag != passiveActiveReal[stmtPos]) {
          // primal return value is currently not updated here (would be the same)
          HandleFactory::template callForwardHandle<PrimalValueTape<TapeTypes> >(statements[stmtPos], 1.0, lhsAdj, passiveActiveReal[stmtPos], indexPos, indices, passivePos, passives, constantPos, constants, primals, adjointData);

#if CODI_EnableVariableAdjointInterfaceInPrimalTapes
          adjointData->setLhsTangent(adjPos); /* Resets the lhs tangent, too */
#else
          adjointData[adjPos] = lhsAdj;
#endif
        }
        stmtPos += 1;
      }
    }

    CODI_INLINE void evaluateStackPrimal(const size_t& startAdjPos, const size_t& endAdjPos,
                                          size_t& constantPos, const size_t& endConstPos, PassiveReal* &constants,
                                          size_t& passivePos, const size_t& endPassivePos, Real* &passives,
                                          size_t& indexPos, const size_t& endIndexPos, Index* &indices,
                                          size_t& stmtPos, const size_t& endStmtPos, Handle* &statements,
                                          StatementInt* &passiveActiveReal) {
      CODI_UNUSED(endConstPos);
      CODI_UNUSED(endPassivePos);
      CODI_UNUSED(endIndexPos);
      CODI_UNUSED(endStmtPos);

      size_t adjPos = startAdjPos;

      while(adjPos < endAdjPos) {
        adjPos += 1;

        if(StatementIntInputTag != passiveActiveReal[stmtPos]) {
          primals[adjPos] = HandleFactory::template callPrimalHandle<PrimalValueTape<TapeTypes> >(statements[stmtPos], passiveActiveReal[stmtPos], indexPos, indices, passivePos, passives, constantPos, constants, primals);
        }

        stmtPos += 1;
      }
    }

  public:

    /**
     * @brief Special evaluation function for the preaccumulation of a tape part.
     *
     * The function just evaluates the tape and does not store the data for the preaccumulation.
     * This function can be used by the tape implementation to reset its state in a more efficient way
     * then it could be programmed from the outside.
     *
     * It has to hold start >= end.
     *
     * @param[in] start The starting position for the reverse evaluation.
     * @param[in]   end The ending position for the reverse evaluation.
     */
    CODI_INLINE void evaluatePreacc(const Position& start, const Position& end) {

      evaluate(start, end);
    }


    /**
     * @brief Special forward evaluation function for the preaccumulation of a tape part.
     *
     * The function just evaluates the tape and does not store the data for the preaccumulation.
     * This function can be used by the tape implementation to reset its state in a more efficient way
     * then it could be programmed from the outside.
     *
     * It has to hold start <= end.
     *
     * @param[in] start The starting position for the reverse evaluation.
     * @param[in]   end The ending position for the reverse evaluation.
     */
    CODI_INLINE void evaluateForwardPreacc(const Position& start, const Position& end) {

      evaluateForward(start, end);
    }

    /**
     * @brief Register a variable as an active variable.
     *
     * The index of the variable is set to a non zero number.
     * @param[in,out] value The value which will be marked as an active variable.
     */
    CODI_INLINE void registerInput(ActiveReal<PrimalValueTape<TapeTypes> >& value) {
      if(isActive()) {
        pushStmtData(value.getGradientData(), value.getValue(), HandleFactory::template createHandle<InputExpr<Real>, PrimalValueTape<TapeTypes> >(), StatementInt(StatementIntInputTag));
      }
    }

    /**
     * @brief Modify the output of an external function such that the tape sees it as an active variable.
     *
     * @param[in,out] value  The output value of the external function.
     * @return Zero
     */
    CODI_INLINE Real registerExtFunctionOutput(ActiveReal<PrimalValueTape<TapeTypes> >& value) {
      registerInput(value);

      return Real();
    }

    /**
     * @brief It is ensured that each output variable has a unique index.
     *
     * @param[in] value  The value will have an unique index that is used by no other variable.
     */
    CODI_INLINE void registerOutput(ActiveReal<PrimalValueTape<TapeTypes> >& value) {
      if(isActive() && value.getGradientData() != 0) {
        Index rhsIndex = value.getGradientData();

        pushCopyHandle(value.getValue(), value.getGradientData(), rhsIndex);
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

  #include "modules/primalValueStaticModule.tpp"
  #undef TAPE_NAME

}
