/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2019 Chair for Scientific Computing (SciComp), TU Kaiserslautern
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
#include "modules/externalFunctionsModule.hpp"
#include "modules/ioModule.hpp"
#include "modules/primalValueModule.hpp"
#include "modules/tapeBaseModule.hpp"
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

    /** @brief The gradient data is just the index type. */
    typedef Index GradientData;

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
  class PrimalValueTape :
      public TapeBaseModule<TapeTypes, PrimalValueTape<TapeTypes>>,
      public PrimalValueModule<TapeTypes, PrimalValueTape<TapeTypes>>,
      public ExternalFunctionModule<TapeTypes, PrimalValueTape<TapeTypes>>,
      public IOModule<TapeTypes, PrimalValueTape<TapeTypes>>,
      public virtual ReverseTapeInterface<typename TapeTypes::Real, typename TapeTypes::Index, typename TapeTypes::GradientValue, PrimalValueTape<TapeTypes>, typename TapeTypes::Position >
  {
  public:

    friend TapeBaseModule<TapeTypes, PrimalValueTape>;  /**< No doc */
    friend PrimalValueModule<TapeTypes, PrimalValueTape>;  /**< No doc */
    friend ::codi::IOModule<TapeTypes, PrimalValueTape>;  /**< No doc */

    CODI_INLINE_REVERSE_TAPE_TYPES(TapeTypes::BaseTypes)

    /** @brief The tape type structure, that defines the basic types. */
    typedef typename TapeTypes::BaseTypes BaseTypes;

    /** @brief The generator for all handles to expression templates. */
    typedef typename TapeTypes::HandleFactory HandleFactory;

    /** @brief The actual type of the stored handles. */
    typedef typename HandleFactory::Handle Handle;

    /** @brief The gradient data is just the index type. */
    typedef typename TapeTypes::GradientData GradientData;

    /** @brief The global position for the tape */
    typedef typename TapeTypes::Position Position;

    /** @brief Vector type of the adjoint vector. Defaults to the template argument but can be switch to a gneral interface. */
    template<typename AdjointData>
    using AdjVecType = typename PrimalValueModule<TapeTypes, PrimalValueTape>::template AdjVecType<AdjointData>;

    /** @brief Adjoint vector interface type. Default vector is the template argument but can be switch to a gneral interface. */
    template<typename AdjointData>
    using AdjVecInterface = typename PrimalValueModule<TapeTypes, PrimalValueTape>::template AdjVecInterface<AdjointData>;

    /** @brief The index handler for the active real's. */
    IndexHandler indexHandler;

    /** @brief Disables code path in CoDiPack that are optimized for Jacobi taping */
    static const bool AllowJacobiOptimization = false;

    /** @brief This tape requires no special primal value handling since the primal value vector is not overwritten. */
    static const bool RequiresPrimalReset = false;

  public:
    /**
     * @brief Creates a tape with the size of zero for the data, statements and external functions.
     */
    PrimalValueTape() :
        TapeBaseModule<TapeTypes, PrimalValueTape<TapeTypes> > (),
        PrimalValueModule<TapeTypes, PrimalValueTape>(),
        ExternalFunctionModule<TapeTypes, PrimalValueTape<TapeTypes> > (),
        IOModule<TapeTypes, PrimalValueTape<TapeTypes> > (),
        indexHandler(MaxStatementIntSize - 1) {
      this->initPrimalValueModule(&indexHandler);
      this->initExtFuncModule(&this->constantValueVector);
      this->initIOModule();
      this->initTapeBaseModule();
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
      this->swapTapeBaseModule(other);
      this->swapPrimalValueModule(other);

      this->extFuncVector.swap(other.extFuncVector);
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

      Index startPos = min(end.inner.inner.inner.inner.inner, this->adjointsSize - 1);
      Index endPos = min(start.inner.inner.inner.inner.inner, this->adjointsSize - 1);

      for(Index i = startPos + 1; i <= endPos; ++i) {
        this->adjoints[i] = GradientValue();
      }
    }
    using TapeBaseModule<TapeTypes, PrimalValueTape>::clearAdjoints;

    /**
     * @brief Set the size of the index and statement data and the primal vector.
     * @param[in] dataSize  The new size of the index vector.
     * @param[in] stmtSize  The new size of the statement vector.
     */
    void resize(const size_t& dataSize, const size_t& stmtSize) {
      this->indexVector.resize(dataSize);
      this->stmtVector.resize(stmtSize);

      this->resizePrimals(stmtSize + 1);
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
      this->stmtVector.reserveItems(1);
      this->stmtVector.setDataAndMove(handle, passiveVariableNumber);
      indexHandler.assignIndex(lhsIndex);

      this->checkPrimalsSize();
      this->primals[lhsIndex] = rhsValue;
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
      ENABLE_CHECK(OptTapeActivity, this->active){
        lhsIndex = rhs.getGradientData();
      } else {
        indexHandler.freeIndex(lhsIndex);
      }
      lhsValue = rhs.getValue();
    }
    using PrimalValueModule<TapeTypes, PrimalValueTape>::store;

  private:

    /**
     * @brief Get the root vector for general data operations.
     *
     * @return The root vector for general data operations.
     */
    CODI_INLINE typename TapeTypes::ExternalFunctionVector& getRootVector() {
      return this->extFuncVector;
    }

    /**
     * @brief Get the root vector for general data operations.
     *
     * @return The root vector for general data operations.
     */
    CODI_INLINE const typename TapeTypes::ExternalFunctionVector& getRootVector() const {
      return this->extFuncVector;
    }

    /**
     * @brief Reset the tape structure to the given position.
     *
     * The state of the tape is then such that all recorded data after this position is
     * no longer used in the evaluation.
     *
     * The allocated memory is not freed. It is used for the next recording.
     *
     * @param[in] pos  The position
     */
    CODI_INLINE void resetInternal(const Position& pos) {
      this->resetExtFunc(pos);
    }

    /**
     * @brief Evaluate the stack from the start to to the end position.
     *
     * It has to hold start >= end.
     *
     * @param[in]       startAdjPos  The starting position for the adjoint evaluation.
     * @param[in]         endAdjPos  The ending position for the adjoint evaluation.
     * @param[in,out]    primalData  The vector of the primal variables.
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
    static CODI_INLINE void evaluateStackReverse(const size_t& startAdjPos, const size_t& endAdjPos, Real* primalData, AdjointData* adjointData,
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
          HandleFactory::template callHandle<PrimalValueTape<TapeTypes> >(statements[stmtPos], 1.0, passiveActiveReal[stmtPos], indexPos, indices, passivePos, passives, constantPos, constants, primalData, adjointData);
#else
          HandleFactory::template callHandle<PrimalValueTape<TapeTypes> >(statements[stmtPos], adj, passiveActiveReal[stmtPos], indexPos, indices, passivePos, passives, constantPos, constants, primalData, adjointData);
#endif
        }
      }
    }

    WRAP_FUNCTION_TEMPLATE(Wrap_evaluateStackReverse, evaluateStackReverse);

    /**
     * @brief Evaluate the stack in the forward mode from the start to to the end position.
     *
     * It has to hold start <= end.
     *
     * @param[in]       startAdjPos  The starting position for the adjoint evaluation.
     * @param[in]         endAdjPos  The ending position for the adjoint evaluation.
     * @param[in,out]    primalData  The vector of the primal variables.
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
    static CODI_INLINE void evaluateStackForward(const size_t& startAdjPos, const size_t& endAdjPos, Real* primalData, AdjointData* adjointData,
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
          HandleFactory::template callForwardHandle<PrimalValueTape<TapeTypes> >(statements[stmtPos], 1.0, lhsAdj, passiveActiveReal[stmtPos], indexPos, indices, passivePos, passives, constantPos, constants, primalData, adjointData);

#if CODI_EnableVariableAdjointInterfaceInPrimalTapes
          adjointData->setLhsTangent(adjPos); /* Resets the lhs tangent, too */
#else
          adjointData[adjPos] = lhsAdj;
#endif
        }
        stmtPos += 1;
      }
    }

    WRAP_FUNCTION_TEMPLATE(Wrap_evaluateStackForward, evaluateStackForward);

    /**
     * @brief Evaluate the stack from the start to to the end position for the primal evaluation.
     *
     * It has to hold start <= end.
     *
     * @param[in]       startAdjPos  The starting position for the adjoint evaluation.
     * @param[in]         endAdjPos  The ending position for the adjoint evaluation.
     * @param[in,out]    primalData  The vector of the primal variables.
     * @param[in,out]   constantPos  The current position in the constant data vector. It will decremented in the method.
     * @param[in]    endConstantPos  The ending position in the constant data vector.
     * @param[in]         constants  The constant values in the rhs expressions.
     * @param[in,out]      indexPos  The current position for the index data. It will decremented in the method.
     * @param[in]       endIndexPos  The ending position for the index data.
     * @param[in]           indices  The indices for the arguments of the rhs.
     * @param[in,out]       stmtPos  The current position in the statement data. It will decremented in the method.
     * @param[in]        endStmtPos  The ending position for statement data.
     * @param[in]        lhsIndices  The indices from the lhs of each statement.
     * @param[in]     storedPrimals  The overwritten primal from the primal vector.
     * @param[in]        statements  The vector with the handles for each statement.
     * @param[in] passiveActiveReal  The number passive values for each statement.
     */
    static CODI_INLINE void evaluateStackPrimal(const size_t& startAdjPos, const size_t& endAdjPos, Real* primalData,
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
          primalData[adjPos] = HandleFactory::template callPrimalHandle<PrimalValueTape<TapeTypes> >(statements[stmtPos], passiveActiveReal[stmtPos], indexPos, indices, passivePos, passives, constantPos, constants, primalData);
        }

        stmtPos += 1;
      }
    }

    WRAP_FUNCTION(Wrap_evaluateStackPrimal, evaluateStackPrimal);


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
    CODI_INLINE void evaluateInternal(const Position& start, const Position& end, AdjointData* adjointData) {

      AdjVecInterface<AdjointData> interface(adjointData, this->primals);
      AdjVecType<AdjointData>* adjVec = this->wrapAdjointVector(interface, adjointData);

      Wrap_evaluateStackReverse<AdjVecType<AdjointData>> evalFunc{};
      auto reverseFunc = &TapeTypes::ConstantValueVector::template evaluateReverse<decltype(evalFunc), Real*&, AdjVecType<AdjointData>*&>;
      this->evaluateExtFunc(start, end, reverseFunc, this->constantValueVector, &interface, evalFunc, this->primals, adjVec);
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
    CODI_INLINE void evaluateForwardInternal(const Position& start, const Position& end, AdjointData* adjointData) {

      AdjVecInterface<AdjointData> interface(adjointData, this->primals);
      AdjVecType<AdjointData>* adjVec = this->wrapAdjointVector(interface, adjointData);

      Wrap_evaluateStackForward<AdjVecType<AdjointData>> evalFunc{};
      auto forwardFunc = &TapeTypes::ConstantValueVector::template evaluateForward<decltype(evalFunc), Real*&, AdjVecType<AdjointData>*&>;
      this->evaluateExtFuncForward(start, end, forwardFunc, this->constantValueVector, &interface, evalFunc, this->primals, adjVec);
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
    CODI_INLINE void evaluatePrimalInternal(const Position& start, const Position& end) {

      AdjVecInterface<GradientValue> interface(this->adjoints, this->primals);

      Wrap_evaluateStackPrimal evalFunc{};
      auto primalFunc = &TapeTypes::ConstantValueVector::template evaluateForward<decltype(evalFunc), Real*&>;
      this->evaluateExtFuncPrimal(start, end, primalFunc, this->constantValueVector, &interface, evalFunc, this->primals);
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

      this->evaluate(start, end);
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

      this->evaluateForward(start, end);
    }

    /**
     * @brief Register a variable as an active variable.
     *
     * The index of the variable is set to a non zero number.
     * @param[in,out] value The value which will be marked as an active variable.
     */
    CODI_INLINE void registerInput(ActiveReal<PrimalValueTape<TapeTypes> >& value) {
      if(this->isActive()) {
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
      if(this->isActive() && value.getGradientData() != 0) {
        Index rhsIndex = value.getGradientData();

        this->pushCopyHandle(value.getValue(), value.getGradientData(), rhsIndex);
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

      this->addTapeBaseValues(values);
      this->addPrimalValueValues(values);
      this->addExtFuncValues(values);

      return values;
    }
  };
}
