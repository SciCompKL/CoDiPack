/**
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015 Chair for Scientific Computing, TU Kaiserslautern
 *
 * This file is part of CoDiPack.
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
 * Authors: Max Sagebaum, Tim Albring, TU Kaiserslautern.
 */

#pragma once

#include <cstddef>
#include <tuple>

#include "../activeReal.hpp"
#include "chunk.hpp"
#include "reverseTapeInterface.hpp"

namespace codi {

  /**
   * @brief Position for the simple tape.
   */
  struct SimpleTapePosition {
    /** @brief The current statement recorded on the tape. */
    size_t stmt;
    /** @brief The current jacobi data recorded on the tape. */
    size_t data;
    /** @brief The current external function recorded on the tape. */
    size_t extFunc;

    /**
     * @brief Simple constructor for convenience.
     * @param[in]    stmt  The current statement recorded on the tape.
     * @param[in]    data  The current jacobi recorded on the tape.
     * @param[in] extFunc  The current external function recorded on the tape.
     */
    SimpleTapePosition(const size_t& stmt, const size_t& data, const size_t& extFunc) :
      stmt(stmt),
      data(data),
      extFunc(extFunc) {}
  };

  /**
   * @brief A tape with a simple implementation and no bounds checking.
   *
   * The SimpleTape implements a fully featured ReverseTapeInterface in a
   * simple fashion. This tape is not intended for simple usage. Actually the
   * tape has no bounds checking, therefore it can produce segmentation faults
   * if it is not used with care.
   *
   * The size of the tape can be set with the resize function and the setExternalFunctionChunkSize.
   *
   * For details on how this tape works please read the general documentation //TODO: Add reference to chapter.
   *
   * Assertions are placed in all the functions such that during development no
   * bounds are overwritten.
   *
   * @tparam      Real  The floating point type used in the ActiveReal.
   * @tparam IndexType  The type for the indexing of the adjoint variables.
   */
  template <typename Real, typename IndexType>
  class SimpleTape : public ReverseTapeInterface<Real, IndexType, SimpleTape<Real, IndexType>, SimpleTapePosition > {
  public:

    /**
     * @brief The type used to store the position of the tape.
     */
    typedef SimpleTapePosition Position;

  private:
    /**
     * @brief The jacobi and index data for the reverse evaluation.
     */
    Chunk2<Real, IndexType> data;
    /**
     * @brief The number of active variables in each statement.
     */
    Chunk1<StatementInt> statements;
    /**
     * @brief The external function data and the position where the external function has been inserted.
     */
    Chunk2<ExternalFunction, Position> externalFunctions;

    /**
     * @brief The adjoint vector.
     */
    Chunk1<Real> adjoints;

    /**
     * @brief Determines if statements are recorded or ignored.
     */
    bool active;

  public:
    /**
     * @brief Creates a tape with the size of zero for the data, statements and external functions.
     */
    SimpleTape() :
      data(0),
      statements(0),
      externalFunctions(0),
      adjoints(1),
      active(false){}

    /**
     * @brief Set the size for the external functions.
     *
     * The method is called this way in order to be compatible with the ChunkTape. It sets the
     * total size of the external functions.
     *
     * @param[in] extChunkSize The new size of
     */
    void setExternalFunctionChunkSize(const size_t& extChunkSize) {
      externalFunctions.resize(extChunkSize);
    }

    /**
     * @brief Set the size of the jacobi and statement data and the adjoint vector.
     * @param[in] dataSize  The new size of the jacobi vector.
     * @param[in] stmtSize  The new size of the statement vector.
     */
    void resize(const size_t& dataSize, const size_t& stmtSize) {
      data.resize(dataSize);
      statements.resize(stmtSize);
      adjoints.resize(stmtSize + 1);
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
      Real gradient; /* This value will not be used */

      ENABLE_CHECK(OptTapeActivity, active){
        assert(ExpressionTraits<Rhs>::maxActiveVariables < data.getUnusedSize());
        /* first store the size of the current stack position and evaluate the
         rhs expression. If there was an active variable on the rhs, update
         the index of the lhs */
        size_t startSize = data.getUsedSize();
        rhs.calcGradient(gradient);
        size_t activeVariables = data.getUsedSize() - startSize;
        if(0 == activeVariables) {
          lhsIndex = 0;
        } else {
          assert(statements.getUsedSize() < statements.size);
          statements.setDataAndMove(std::make_tuple((StatementInt)activeVariables));
          lhsIndex = statements.getUsedSize();
        }
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
    inline void store(Real& lhsValue, IndexType& lhsIndex, const ActiveReal<Real, SimpleTape<Real, IndexType> >& rhs) {
      ENABLE_CHECK(OptTapeActivity, active){
        lhsIndex = rhs.getGradientData();
      } else {
        lhsIndex = 0;
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
      lhsIndex = 0;
      lhsValue = rhs;
    }

    /**
     * @brief Stores the jacobi with the value 1.0 on the tape if the index is active.
     *
     * @param[in] gradient Not used in this implementation.
     * @param[in]    value Not used in this implementation.
     * @param[in]    index Used to check if the variable is active.
     */
    inline void pushJacobi(Real& gradient, const Real& value, const IndexType& index) {
      CODI_UNUSED(gradient);
      CODI_UNUSED(value);

      if(0 != index) {
        assert(data.getUsedSize() < data.size);

        data.setDataAndMove(std::make_tuple(1.0, index));
      }
    }

    /**
     * @brief Stores the jacobi on the tape if the index is active.
     *
     * @param[in] gradient Not used in this implementation.
     * @param[in]   jacobi Stored on the tape if the variable is active.
     * @param[in]    value Not used in this implementation.
     * @param[in]    index Used to check if the variable is active.
     */
    inline void pushJacobi(Real& gradient, const Real& jacobi, const Real& value, const IndexType& index) {
      CODI_UNUSED(gradient);
      CODI_UNUSED(value);

      if(0 != index) {
        ENABLE_CHECK(OptIgnoreInvalidJacobies, isfinite(jacobi)) {
          ENABLE_CHECK(OptJacobiIsZero, 0.0 != jacobi) {
            assert(data.getUsedSize() < data.size);

            data.setDataAndMove(std::make_tuple(jacobi, index));
          }
        }
      }
    }

    /**
     * @brief Set the index to zero.
     * @param[in] value Not used in this implementation.
     * @param[out] index The index of the active type.
     */
    inline void initGradientData(Real& value, IndexType& index) {
      CODI_UNUSED(value);
      index = 0;
    }

    /**
     * @brief Does nothing.
     * @param[in] value Not used in this implementation.
     * @param[in] index Not used in this implementation.
     */
    inline void destroyGradientData(Real& value, IndexType& index) {
      CODI_UNUSED(value);
      CODI_UNUSED(index);
      /* nothing to do */
    }


    /**
     * @brief Set the gradient value of the corresponding index.
     *
     * If the index 0 is the inactive indicator and is ignored.
     *
     * @param[in]    index  The index of the active type.
     * @param[in] gradient  The new value for the gradient.
     */
    void setGradient(IndexType& index, const Real& gradient) {
      if(0 != index) {
        this->gradient(index) = gradient;
      }
    }

    /**
     * @brief Get the gradient value of the corresponding index.
     *
     * @param[in] index The index of the active type.
     * @return The gradient value corresponding to the given index.
     */
    inline Real getGradient(const IndexType& index) const {
      assert((size_t)index < statements.size);
      return adjoints.data[index];
    }

    /**
     * @brief Get a reference to the gradient value of the corresponding index.
     *
     * An index of 0 will raise an assert exception.
     *
     * @param[in] index The index of the active type.
     * @return The reference to the gradient data.
     */
    inline Real& gradient(IndexType& index) {
      assert((size_t)index < statements.size);
      assert(0 != index);

      return adjoints.data[index];
    }

    /**
     * @brief Get the current position of the tape.
     *
     * The position can be used to reset the tape to that position or to
     * evaluate only parts of the tape.
     * @return The current position of the tape.
     */
    inline Position getPosition() {
      return Position(statements.getUsedSize(), data.getUsedSize(), externalFunctions.getUsedSize());
    }

    /**
     * @brief Reset the tape to the given position.
     *
     * @param[in] pos Reset the state of the tape to the given position.
     */
    inline void reset(const Position& pos) {
      assert(pos.stmt < statements.size);
      assert(pos.data < data.size);
      assert(pos.extFunc < externalFunctions.size);

      for(size_t i = pos.stmt; i <= statements.getUsedSize(); ++i) {
        adjoints.data[i] = 0.0;
      }

      for(size_t i = pos.extFunc; i < externalFunctions.getUsedSize(); ++i) {
        externalFunctions.data1[i].deleteData();
      }

      statements.setUsedSize(pos.stmt);
      data.setUsedSize(pos.data);
      externalFunctions.setUsedSize(pos.extFunc);
    }

    /**
     * @brief Reset the tape to its initial state.
     */
    inline void reset() {
      reset(Position(0,0,0));
    }

    /**
     * @brief Sets all adjoint/gradients to zero.
     */
    inline void clearAdjoints(){
      for(size_t i = 0; i <= statements.getUsedSize(); ++i) {
        adjoints.data[i] = 0.0;
      }
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
    inline void evaluateStack(const Position& start, const Position& end) {
      Position curPos = start;

      while(curPos.stmt > end.stmt) {
        const Real& adj = adjoints.data[curPos.stmt];
        --curPos.stmt;
        const StatementInt& activeVariables = statements.data[curPos.stmt];
        ENABLE_CHECK(OptZeroAdjoint, adj != 0.0){
          for(StatementInt curVar = 0; curVar < activeVariables; ++curVar) {
            --curPos.data;

            adjoints.data[data.data2[curPos.data]] += adj * data.data1[curPos.data];
          }
        } else {
          curPos.data -= activeVariables;
        }
      }
    }

  public:
    /**
     * @brief Perform the adjoint evaluation from start to end.
     *
     * It has to hold start >= end.
     *
     * @param[in] start  The starting position for the adjoint evaluation.
     * @param[in]   end  The ending position for the adjoint evaluation.
     */
    inline void evaluate(const Position& start, const Position& end) {
      assert(start.data >= end.data);
      assert(start.stmt >= end.stmt);
      assert(start.extFunc >= end.extFunc);

      Position curPos = start;


      for(size_t curExtFunc = start.extFunc; curExtFunc > end.extFunc; /* decrement is done inside the loop */) {
        --curExtFunc; // decrement of loop variable

        ExternalFunction& extFunc = externalFunctions.data1[curExtFunc];
        const Position& extFuncPos = externalFunctions.data2[curExtFunc];

        // always evaluate the stack to the point of the external function
        evaluateStack(curPos, extFuncPos);

        extFunc.evaluate();

        curPos = extFuncPos;
      }

      // Iterate over the reminder also covers the case if the there are no external functions
      evaluateStack(curPos, end);
    }

    /**
     * @brief Perform the adjoint evaluation from the current position to the initial position.
     */
    inline void evaluate() {
      evaluate(getPosition(), Position(0,0,0));
    }

    /**
     * @brief Register a variable as an active variable.
     *
     * The index of the variable is set to a non zero number.
     * @param[inout] value The value which will be marked as an active variable.
     */
    inline void registerInput(ActiveReal<Real, SimpleTape<Real, IndexType> >& value) {
      assert(statements.getUsedSize() < statements.size);

      statements.setDataAndMove(std::make_tuple((StatementInt) 0));
      value.getGradientData() = statements.getUsedSize();
    }

    /**
     * @brief Not needed in this implementation.
     *
     * @param[in] value Not used.
     */
    inline void registerOutput(ActiveReal<Real, SimpleTape<Real, IndexType> >& value) {
      CODI_UNUSED(value);
      /* do nothing */
    }

    /**
     * @brief Start recording.
     */
    inline void setActive(){
      active = true;
    }

    /**
     * @brief Stop recording.
     */
    inline void setPassive(){
      active = false;
    }

    /**
     * @brief Check if the tape is active.
     * @return true if the tape is active.
     */
    inline bool isActive(){
      ENABLE_CHECK(OptTapeActivity, true) {
        // default branch will return the tape activity
        return active;
      } else {
        // if we do not check for the tape activity, the tape is always active
        return true;
      }
    }

    /**
     * @brief Add an external function with a void handle as user data.
     *
     * The data handle provided to the tape is considered in possession of the tape. The tape will now be responsible to
     * free the handle. For this it will use the delete function provided by the user.
     *
     * @param[in] extFunc  The external function which is called by the tape.
     * @param[inout] data  The data for the external function. The tape takes ownership over the data.
     * @param[in] delData  The delete function for the data.
     */
    void pushExternalFunctionHandle(ExternalFunction::CallFunction extFunc, void* data, ExternalFunction::DeleteFunction delData){
      pushExternalFunctionHandle(ExternalFunction(extFunc, data, delData));
    }


    /**
     * @brief Add an external function with a specific data type.
     *
     * The data pointer provided to the tape is considered in possession of the tape. The tape will now be responsible to
     * free the data. For this it will use the delete function provided by the user.
     *
     * @param[in] extFunc  The external function which is called by the tape.
     * @param[inout] data  The data for the external function. The tape takes ownership over the data.
     * @param[in] delData  The delete function for the data.
     */
    template<typename Data>
    void pushExternalFunction(typename ExternalFunctionDataHelper<Data>::CallFunction extFunc, Data* data, typename ExternalFunctionDataHelper<Data>::DeleteFunction delData){
      pushExternalFunctionHandle(ExternalFunctionDataHelper<Data>::createHandle(extFunc, data, delData));\
    }

  private:
    /**
     * @brief Private common method to add to the external function stack.
     *
     * @param[in] function The external function structure to push.
     */
    void pushExternalFunctionHandle(const ExternalFunction& function){
      assert(0 != externalFunctions.getUnusedSize());
      externalFunctions.setDataAndMove(std::make_tuple(function, getPosition()));
    }

  };
}
