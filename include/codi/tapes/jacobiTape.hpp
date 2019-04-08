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
#include <iostream>
#include <iomanip>
#include <map>
#include <tuple>

#include "../activeReal.hpp"
#include "../typeFunctions.hpp"
#include "chunk.hpp"
#include "chunkVector.hpp"
#include "externalFunctions.hpp"
#include "modules/externalFunctionsModule.hpp"
#include "modules/ioModule.hpp"
#include "modules/jacobiModule.hpp"
#include "modules/statementModule.hpp"
#include "modules/tapeBaseModule.hpp"
#include "reverseTapeInterface.hpp"
#include "singleChunkVector.hpp"
#include "../tapeTypes.hpp"
#include "../tools/tapeValues.hpp"

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
   * @tparam        RTT  The basic type defintions for the tape. Need to define everything from ReverseTapeTypes.
   * @tparam DataVector  The data manager for the chunks. Needs to implement a ChunkVector interface.
   */
  template <typename RTT, template<typename, typename> class DataVector>
  struct JacobiTapeTypes {

    CODI_INLINE_REVERSE_TAPE_TYPES(RTT)

    /** @brief The tape type structure, that defines the basic types. */
    typedef RTT BaseTypes;

    /** @brief The data for each statement. */
    typedef Chunk1<StatementInt> StatementChunk;
    /** @brief The chunk vector for the statement data. */
    typedef DataVector<StatementChunk, IndexHandler> StatementVector;

    /** @brief The data for the jacobies of each statement */
    typedef Chunk2< Real, typename IndexHandler::Index> JacobiChunk;
    /** @brief The chunk vector for the jacobi data. */
    typedef DataVector<JacobiChunk, StatementVector> JacobiVector;

    /** @brief The data for the external functions. */
    typedef Chunk2<ExternalFunction,typename JacobiVector::Position> ExternalFunctionChunk;
    /** @brief The chunk vector for the external  function data. */
    typedef DataVector<ExternalFunctionChunk, JacobiVector> ExternalFunctionVector;

    /** @brief The position for all the different data vectors. */
    typedef typename ExternalFunctionVector::Position Position;

    /** @brief The gradient data is just the index type. */
    typedef Index GradientData;

    /** @brief The name of the tape as a string. */
    constexpr static const char* tapeName = "JacobiTape";

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
  class JacobiTape final :
      public TapeBaseModule<TapeTypes, JacobiTape<TapeTypes>>,
      public JacobiModule<TapeTypes, JacobiTape<TapeTypes>>,
      public StatementModule<TapeTypes, JacobiTape<TapeTypes>>,
      public ExternalFunctionModule<TapeTypes, JacobiTape<TapeTypes>>,
      public IOModule<TapeTypes, JacobiTape<TapeTypes>>,
      public virtual ReverseTapeInterface<typename TapeTypes::Real, typename TapeTypes::Index, typename TapeTypes::GradientValue, JacobiTape<TapeTypes>, typename TapeTypes::Position >
  {
  public:

    friend TapeBaseModule<TapeTypes, JacobiTape>;  /**< No doc */
    friend StatementModule<TapeTypes, JacobiTape>;  /**< No doc */
    friend ::codi::IOModule<TapeTypes, JacobiTape>;  /**< No doc */

    CODI_INLINE_REVERSE_TAPE_TYPES(TapeTypes::BaseTypes)

    /** @brief The tape type structure, that defines all types for the tape. */
    typedef typename TapeTypes::BaseTypes BaseTypes;

    /** @brief The gradient data is just the index type. */
    typedef typename TapeTypes::GradientData GradientData;

    /** @brief The global position for the tape */
    typedef typename TapeTypes::Position Position;

    /** @brief The index handler for the active real's. */
    IndexHandler indexHandler;

    /** @brief Enables code path in CoDiPack that are optimized for Jacobi taping */
    static const bool AllowJacobiOptimization = true;

    /** @brief This tape requires no primal value handling. */
    static const bool RequiresPrimalReset = false;

  public:
    /**
     * @brief Creates a tape with the default chunk sizes for the data, statements and
     * external functions defined in the configuration.
     */
    JacobiTape() :
      TapeBaseModule<TapeTypes, JacobiTape<TapeTypes> > (),
      JacobiModule<TapeTypes, JacobiTape>(),
      StatementModule<TapeTypes, JacobiTape<TapeTypes> > (),
      ExternalFunctionModule<TapeTypes, JacobiTape<TapeTypes> > (),
      IOModule<TapeTypes, JacobiTape<TapeTypes> > (),
      indexHandler(0) {
      this->initStmtModule(&indexHandler);
      this->initJacobiModule(&this->stmtVector);
      this->initExtFuncModule(&this->jacobiVector);
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
    void swap(JacobiTape& other) {
      this->swapTapeBaseModule(other);

      this->extFuncVector.swap(other.extFuncVector);
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
      this->resizeJacobi(dataSize);
      this->resizeStmt(statementSize);
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
      Index startPos = min(end.inner.inner.inner, this->adjointsSize - 1);
      Index endPos = min(start.inner.inner.inner, this->adjointsSize - 1);

      for(Index i = startPos + 1; i <= endPos; ++i) {
        this->adjoints[i] = GradientValue();
      }
    }

    using TapeBaseModule<TapeTypes, JacobiTape>::clearAdjoints;

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
    CODI_INLINE void store(Real& lhsValue, Index& lhsIndex, const ActiveReal<JacobiTape<TapeTypes> >& rhs) {
      ENABLE_CHECK (OptTapeActivity, this->active){
        lhsIndex = rhs.getGradientData();
      } else {
        indexHandler.freeIndex(lhsIndex);
      }
      lhsValue = rhs.getValue();
    }

    using StatementModule<TapeTypes, JacobiTape>::store;

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
     * @brief The callback method for the push of statement data.
     *
     * The method is called by the statement module to push the
     * statements on the tape.
     *
     * @param[in] numberOfArguments  The number of arguments in the statements that have been pushed as jacobies.
     * @param[in]          lhsIndex  The index of the lhs value of the operation.
     */
    CODI_INLINE void pushStmtData(const StatementInt& numberOfArguments, const Index& lhsIndex) {
      CODI_UNUSED(lhsIndex);

      this->stmtVector.setDataAndMove(numberOfArguments);
    }

    /**
     * @brief Implementation of the AD stack evaluation.
     *
     * It has to hold startAdjPos >= endAdjPos.
     *
     * @param[in]     startAdjPos  The starting point in the expression evaluation.
     * @param[in]       endAdjPos  The ending point in the expression evaluation.
     * @param[in,out] adjointData  The vector of the adjoint variables.
     * @param[in,out]     dataPos  The current position in the jacobi and index vector. This value is used in the next invocation of this method.
     * @param[in]      endDataPos  The end position in the jacobi and index vector.
     * @param[in]        jacobies  The pointer to the jacobi vector.
     * @param[in]         indices  The pointer to the index vector
     * @param[in,out]     stmtPos  The current position in the statement vector. This value is used in the next invocation of this method.
     * @param[in]      endStmtPos  The end position in the statement vector.
     * @param[in]      statements  The pointer to the statement vector.
     *
     * @tparam AdjointData The data for the adjoint vector it needs to support add, multiply and comparison operations.
     */
    template<typename AdjointData>
    static CODI_INLINE void evaluateStackReverse(const size_t& startAdjPos, const size_t& endAdjPos, AdjointData* adjointData,
                                      size_t& dataPos, const size_t& endDataPos, Real* &jacobies, Index* &indices,
                                      size_t& stmtPos, const size_t& endStmtPos, StatementInt* &statements) {

      CODI_UNUSED(endDataPos);
      CODI_UNUSED(endStmtPos);

      size_t adjPos = startAdjPos;

      while(adjPos > endAdjPos) {
        --stmtPos;

        const AdjointData adj = adjointData[adjPos];
        if(ZeroAdjointReverse && StatementIntInputTag != statements[stmtPos]) {
          adjointData[adjPos] = AdjointData();
        }
        --adjPos;

#if CODI_AdjointHandle_Jacobi_Reverse
        handleReverseEval(adj, adjPos + 1);
#endif

        if(StatementIntInputTag != statements[stmtPos]) {
          JacobiModule<TapeTypes, JacobiTape>::incrementAdjoints(adj, adjointData, statements[stmtPos], dataPos, jacobies, indices);
        }
      }
    }

    WRAP_FUNCTION_TEMPLATE(Wrap_evaluateStackReverse, evaluateStackReverse);

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
    CODI_INLINE void evaluateInternal(const Position& start, const Position& end, AdjointData* adjointData) {

      Wrap_evaluateStackReverse<AdjointData> evalFunc{};
      auto reverseFunc = &TapeTypes::JacobiVector::template evaluateReverse<decltype(evalFunc), AdjointData*&>;

      AdjointInterfaceImpl<Real, Index, AdjointData> interface(adjointData);

      this->evaluateExtFunc(start, end, reverseFunc, this->jacobiVector, &interface, evalFunc, adjointData);
    }

    /**
     * @brief Implementation of the AD stack evaluation.
     *
     * It has to hold startAdjPos <= endAdjPos.
     *
     * @param[in]     startAdjPos  The starting point in the expression evaluation.
     * @param[in]       endAdjPos  The ending point in the expression evaluation.
     * @param[in,out] adjointData  The vector of the adjoint variables.
     * @param[in,out]     dataPos  The current position in the jacobi and index vector. This value is used in the next invocation of this method.
     * @param[in]      endDataPos  The end position in the jacobi and index vector.
     * @param[in]        jacobies  The pointer to the jacobi vector.
     * @param[in]         indices  The pointer to the index vector
     * @param[in,out]     stmtPos  The current position in the statement vector. This value is used in the next invocation of this method.
     * @param[in]      endStmtPos  The end position in the statement vector.
     * @param[in]      statements  The pointer to the statement vector.
     *
     * @tparam AdjointData The data for the adjoint vector it needs to support add, multiply and comparison operations.
     */
    template<typename AdjointData>
    static CODI_INLINE void evaluateStackForward(const size_t& startAdjPos, const size_t& endAdjPos, AdjointData* adjointData,
                                          size_t& dataPos, const size_t& endDataPos, Real* &jacobies, Index* &indices,
                                          size_t& stmtPos, const size_t& endStmtPos, StatementInt* &statements) {
      CODI_UNUSED(endDataPos);
      CODI_UNUSED(endStmtPos);

      size_t adjPos = startAdjPos;

      while(adjPos < endAdjPos) {
        ++adjPos;
        AdjointData adj = AdjointData();

        if(StatementIntInputTag != statements[stmtPos]) {
          JacobiModule<TapeTypes, JacobiTape>::incrementTangents(adj, adjointData, statements[stmtPos], dataPos, jacobies, indices);
          adjointData[adjPos] = adj;
        }

        ++stmtPos;
      }
    }

    WRAP_FUNCTION_TEMPLATE(Wrap_evaluateStackForward, evaluateStackForward);

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
    CODI_INLINE void evaluateForwardInternal(const Position& start, const Position& end, AdjointData* adjointData) {

      Wrap_evaluateStackForward<AdjointData> evalFunc{};
      auto forwardFunc = &TapeTypes::JacobiVector::template evaluateForward<decltype(evalFunc), AdjointData*&>;

      AdjointInterfaceImpl<Real, Index, AdjointData> interface(adjointData);

      this->evaluateExtFuncForward(start, end, forwardFunc, this->jacobiVector, &interface, evalFunc, adjointData);
    }

    /**
     * @brief Internal function for input variable registration.
     *
     * The index of the variable is set to a non zero number.
     *
     * @param[inout] value  Unused
     * @param[out]   index  Is assigned a non zero value.
     */
    CODI_INLINE void registerInputInternal(Real& value, Index& index) {
      CODI_UNUSED(value);

      this->stmtVector.reserveItems(1);
      this->stmtVector.setDataAndMove((StatementInt)StatementIntInputTag);

      index = indexHandler.createIndex();
    }

    /**
     * @brief Internal function for output variable registration.
     *
     * The index of the variable is set to a non zero number.
     *
     * @param[inout] value  Unused
     * @param[out]   index  Is assigned a non zero value.
     */
    CODI_INLINE void registerOutputInternal(const Real& value, Index& index) {
      CODI_UNUSED(value);

      ENABLE_CHECK(OptCheckZeroIndex, 0 != index) {
        this->stmtVector.reserveItems(1);
        this->jacobiVector.reserveItems(1);

        this->stmtVector.setDataAndMove((StatementInt)1);
        this->jacobiVector.setDataAndMove(1.0, index);

        index = indexHandler.createIndex();
      }
    }

  public:
    /**
     * @brief Register a variable as an active variable.
     *
     * The index of the variable is set to a non zero number.
     * @param[in,out] value The value which will be marked as an active variable.
     */
    CODI_INLINE void registerInput(ActiveReal<JacobiTape<TapeTypes> >& value) {
      registerInputInternal(value.value(), value.getGradientData());
    }

    /**
     * @brief Modify the output of an external function such that the tape sees it as an active variable.
     *
     * @param[in,out] value  The output value of the external function.
     * @return Zero
     */
    CODI_INLINE Real registerExtFunctionOutput(ActiveReal<JacobiTape<TapeTypes> >& value) {
      registerInput(value);

      return Real();
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
      registerOutputInternal(value.getValue(), value.getGradientData());
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
      this->addStmtValues(values);
      this->addJacobiValues(values);
      this->addExtFuncValues(values);
      indexHandler.addValues(values);

      return values;
    }
  };
}
