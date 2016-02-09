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

#include <cstddef>
#include <iomanip>
#include <tuple>

#include "../activeReal.hpp"
#include "chunk.hpp"
#include "indices/reuseIndexHandler.hpp"
#include "reverseTapeInterface.hpp"

/**
 * @brief Global namespace for CoDiPack - Code Differentiation Package
 */
namespace codi {

  /**
   * @brief Position for the simple tape.
   */
  struct SimpleIndexTapePosition {
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
    SimpleIndexTapePosition(const size_t& stmt, const size_t& data, const size_t& extFunc) :
      stmt(stmt),
      data(data),
      extFunc(extFunc) {}
  };

  /**
   * @brief A tape with a simple implementation and no bounds checking.
   *
   * The SimpleIndexTape implements a fully featured ReverseTapeInterface in a
   * simple fashion. This tape is not intended for simple usage. Actually the
   * tape has no bounds checking, therefore it can produce segmentation faults
   * if it is not used with care.
   *
   * The size of the tape can be set with the resize function and the setExternalFunctionChunkSize.
   *
   * For details on how this tape works please read the general documentation //TODO: Add reference to chapter.
   *
   * The tape also uses the index manager ReuseIndexHandler to reuse the indices that are deleted.
   * That means that ActiveReal's which use this tape need to be copied by usual means and deleted after
   * the are no longer used. No c-like memory operations like memset and memcpy should be applied
   * to these types.
   *
   * Assertions are placed in all the functions such that during development no
   * bounds are overwritten.
   *
   * @tparam      Real  The floating point type used in the ActiveReal.
   * @tparam IndexType  The type for the indexing of the adjoint variables.
   */
  template <typename Real, typename IndexType>
  class SimpleIndexTape : public ReverseTapeInterface<Real, IndexType, SimpleIndexTape<Real, IndexType>, SimpleIndexTapePosition > {
  public:

    /**
     * @brief The type used to store the position of the tape.
     */
    typedef SimpleIndexTapePosition Position;

  private:
    /**
     * @brief The jacobi and index data for the reverse evaluation.
     */
    Chunk2<Real, IndexType> data;
    /**
     * @brief The number of active variables in each statement and the index on the lhs.
     */
    Chunk2<StatementInt, IndexType> statements;
    /**
     * @brief The external function data and the position where the external function has been inserted.
     */
    Chunk2<ExternalFunction, Position> externalFunctions;

    /**
     * @brief The adjoint vector.
     */
    Chunk1<Real> adjoints;

    ReuseIndexHandler<IndexType> indexHandler;

    /**
     * @brief Determines if statements are recorded or ignored.
     */
    bool active;

  public:
    /**
     * @brief Creates a tape with the size of zero for the data, statements and external functions.
     */
    SimpleIndexTape() :
      data(0),
      statements(0),
      externalFunctions(0),
      adjoints(1),
      indexHandler(),
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
     * @brief Set the size of the adjoint vector.
     *
     * @param[in] adjointsSize The new size for the adjoint vector.
     */
    void setAdjointsSize(const size_t& adjointsSize) {
      adjoints.resize(adjointsSize);
    }

    /**
     * @brief Return the number of used statements.
     * @return The number of used statements.
     */
    size_t getUsedStatementsSize() {
      return statements.getUsedSize();
    }

    /**
     * @brief Return the number of used data entries.
     * @return The number of used data entries.
     */
    size_t getUsedDataEntriesSize() {
      return data.getUsedSize();
    }

    /**
     * @brief Get the current size of the adjoint vector.
     *
     * @return The size of the current adjoint vector.
     */
    size_t getAdjointsSize() {
      return indexHandler.getMaximumGlobalIndex() + 1;
    }

    /**
     * @brief Set the size of the jacobi and statement data and the adjoint vector.
     * @param[in] dataSize  The new size of the jacobi vector.
     * @param[in] stmtSize  The new size of the statement vector.
     */
    void resize(const size_t& dataSize, const size_t& stmtSize) {
      data.resize(dataSize);
      statements.resize(stmtSize);
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
      void* null = NULL;
      ENABLE_CHECK(OptTapeActivity, active){
        assert(ExpressionTraits<Rhs>::maxActiveVariables < data.getUnusedSize());
        /* first store the size of the current stack position and evaluate the
         rhs expression. If there was an active variable on the rhs, update
         the index of the lhs */
        size_t startSize = data.getUsedSize();
        rhs.template calcGradient<void*>(null);
        size_t activeVariables = data.getUsedSize() - startSize;
        ENABLE_CHECK(OptCheckEmptyStatements, 0 != activeVariables) {
          indexHandler.checkIndex(lhsIndex);
          assert(lhsIndex < (IndexType)adjoints.size);

          assert(statements.getUsedSize() < statements.size);
          statements.setDataAndMove(std::make_tuple((StatementInt)activeVariables, lhsIndex));
        } else {
          indexHandler.freeIndex(lhsIndex);
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
    inline void store(Real& lhsValue, IndexType& lhsIndex, const ActiveReal<Real, SimpleIndexTape<Real, IndexType> >& rhs) {
      ENABLE_CHECK(OptTapeActivity, active){
        ENABLE_CHECK(OptCheckZeroIndex, 0 != rhs.getGradientData()) {
          indexHandler.checkIndex(lhsIndex);
          assert(lhsIndex < (IndexType)adjoints.size);

          assert(statements.getUsedSize() < statements.size);
          assert(1 <= data.getUnusedSize());
          this->data.setDataAndMove(std::make_tuple(1.0, rhs.getGradientData()));
          statements.setDataAndMove(std::make_tuple((StatementInt)1, lhsIndex));
        } else {
          indexHandler.freeIndex(lhsIndex);
        }
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
      indexHandler.freeIndex(lhsIndex);
      lhsValue = rhs;
    }

    /**
     * @brief Manual store routine.
     *
     * Use this routine to add a statement if the corresponding jacobi entries will be manually pushed onto the tape.
     *
     * The Jacobi entries must be pushed immediately after calling this routine using pushJacobi.
     *
     * @param[out]   lhsIndex    The gradient data of the lhs.
     * @param[in]        size    The number of Jacobi entries.
     */
    inline void store(IndexType& lhsIndex, StatementInt size) {
      ENABLE_CHECK (OptTapeActivity, active){
        assert(size < data.getUnusedSize());
        indexHandler.checkIndex(lhsIndex);
        assert(lhsIndex < (IndexType)adjoints.size);
        assert(statements.getUsedSize() < statements.size);
        statements.setDataAndMove(std::make_tuple(size, lhsIndex));
      }
    }

    /**
     * @brief Stores the jacobi with the value 1.0 on the tape if the index is active.
     *
     * @param[in]     data Not used in this implementation.
     * @param[in]    value Not used in this implementation.
     * @param[in]    index Used to check if the variable is active.
     *
     * @tparam Data  The type of the data for the tape.
     */
    template<typename Data>
    inline void pushJacobi(Data& data, const Real& value, const IndexType& index) {
      CODI_UNUSED(data);
      CODI_UNUSED(value);

      ENABLE_CHECK(OptCheckZeroIndex, 0 != index) {
        assert(this->data.getUsedSize() < this->data.size);

        this->data.setDataAndMove(std::make_tuple(1.0, index));
      }
    }

    /**
     * @brief Stores the jacobi on the tape if the index is active.
     *
     * @param[in]     data Not used in this implementation.
     * @param[in]   jacobi Stored on the tape if the variable is active.
     * @param[in]    value Not used in this implementation.
     * @param[in]    index Used to check if the variable is active.
     *
     * @tparam Data  The type of the data for the tape.
     */
    template<typename Data>
    inline void pushJacobi(Data& data, const Real& jacobi, const Real& value, const IndexType& index) {
      CODI_UNUSED(data);
      CODI_UNUSED(value);

      ENABLE_CHECK(OptCheckZeroIndex, 0 != index) {
        ENABLE_CHECK(OptIgnoreInvalidJacobies, isfinite(jacobi)) {
          ENABLE_CHECK(OptJacobiIsZero, 0.0 != jacobi) {
            assert(this->data.getUsedSize() < this->data.size);

            this->data.setDataAndMove(std::make_tuple(jacobi, index));
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
      * @brief Frees the index.
      *
      * @param[in] value Not used in this implementation.
      * @param[in] index The index is given to the index handler.
      */
     inline void destroyGradientData(Real& value, IndexType& index) {
      CODI_UNUSED(value);

      indexHandler.freeIndex(index);
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
      assert((size_t)index < adjoints.size);
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
      assert((size_t)index < adjoints.size);
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

      clearAdjoints();

      for(size_t i = pos.extFunc; i < externalFunctions.getUsedSize(); ++i) {
        externalFunctions.data1[i].deleteData();
      }

      statements.setUsedSize(pos.stmt);
      data.setUsedSize(pos.data);
      externalFunctions.setUsedSize(pos.extFunc);

      indexHandler.reset();
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
      for(IndexType i = 0; i <= indexHandler.getMaximumGlobalIndex(); ++i) {
        adjoints.data[i] = 0.0;
      }
    }

    /**
     * @brief Does nothing because the indices are not connected to the positions.
     *
     * @param[in] start Not used
     * @param[in] end Not used
     */
    inline void clearAdjoints(const Position& start, const Position& end){
      CODI_UNUSED(start);
      CODI_UNUSED(end);
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

      // Usage of pointers increases evaluation speed
      Real* adjoint = adjoints.data;
      StatementInt* argumentCount = statements.data1;
      IndexType* lhsIndex = statements.data2;
      Real* jacobi = data.data1;
      IndexType* rhsIndex = data.data2;
      while(curPos.stmt > end.stmt) {
        --curPos.stmt;

        const IndexType& index = lhsIndex[curPos.stmt];
        const Real adj = adjoint[index];
        adjoints.data[index] = 0.0;
        const StatementInt& activeVariables = argumentCount[curPos.stmt];
        ENABLE_CHECK(OptZeroAdjoint, adj != 0.0){
          for(StatementInt curVar = 0; curVar < activeVariables; ++curVar) {
            --curPos.data;

            adjoint[rhsIndex[curPos.data]] += adj * jacobi[curPos.data];
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
      assert((IndexType)adjoints.size > indexHandler.getMaximumGlobalIndex());

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
    inline void registerInput(ActiveReal<Real, SimpleIndexTape<Real, IndexType> >& value) {
      indexHandler.checkIndex(value.getGradientData());
    }

    /**
     * @brief Not needed in this implementation.
     *
     * @param[in] value Not used.
     */
    inline void registerOutput(ActiveReal<Real, SimpleIndexTape<Real, IndexType> >& value) {
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
      pushExternalFunctionHandle(ExternalFunctionDataHelper<Data>::createHandle(extFunc, data, delData));
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

  public:
    /**
     * @brief Prints statistics about the tape on the screen
     *
     * Prints information such as stored statements/adjoints and memory usage on screen.
     */
    void printStatistics(){
      const double BYTE_TO_MB = 1.0/1024.0/1024.0;

      size_t nAdjoints      = (size_t)indexHandler.getMaximumGlobalIndex() + 1;
      size_t memoryAdjoints = (double)nAdjoints * (double)sizeof(Real) * BYTE_TO_MB;

      size_t totalStmts    = statements.getUsedSize();
      double  memoryUsedStmts = (double)totalStmts*((double)sizeof(StatementInt) + sizeof(IndexType))* BYTE_TO_MB;
      double  memoryAllocStmts= ((double)statements.getUnusedSize()+(double)statements.getUsedSize())
                                *((double)sizeof(StatementInt) + sizeof(IndexType))* BYTE_TO_MB;
      size_t totalData    = data.getUsedSize();
      double  memoryUsedData = (double)totalData*(double)(sizeof(Real)+sizeof(IndexType))* BYTE_TO_MB;
      double  memoryAllocData= ((double)data.getUsedSize()+(double)data.getUnusedSize())
                                *(double)(sizeof(Real)+sizeof(IndexType))* BYTE_TO_MB;

      size_t maximumGlobalIndex     = (size_t)indexHandler.getMaximumGlobalIndex();
      size_t storedIndices          = (size_t)indexHandler.getNumberStoredIndices();
      size_t currentLiveIndices     = (size_t)indexHandler.getCurrentIndex() - indexHandler.getNumberStoredIndices();

      double memoryStoredIndices    = (double)storedIndices*(double)(sizeof(IndexType)) * BYTE_TO_MB;
      double memoryAllocatedIndices = (double)indexHandler.getNumberAllocatedIndices()*(double)(sizeof(IndexType)) * BYTE_TO_MB;

      size_t nExternalFunc = externalFunctions.getUsedSize();

      std::cout << std::endl
                << "---------------------------------------------" << std::endl
                << "CoDi Tape Statistics (SimpleIndexReuseTape)  " << std::endl
                << "---------------------------------------------" << std::endl
                << "Statements " << std::endl
                << "---------------------------------------------" << std::endl
                << "  Total Number:       " << std::setw(10) << totalStmts   << std::endl
                << "  Memory allocated:   " << std::setiosflags(std::ios::fixed)
                                            << std::setprecision(2)
                                            << std::setw(10)
                                            << memoryAllocStmts << " MB" << std::endl
                << "  Memory used:        " << std::setiosflags(std::ios::fixed)
                                            << std::setprecision(2)
                                            << std::setw(10)
                                            << memoryUsedStmts << " MB" << std::endl
                << "---------------------------------------------" << std::endl
                << "Jacobi entries "                       << std::endl
                << "---------------------------------------------" << std::endl
                << "  Total Number:       " << std::setw(10) << totalData   << std::endl
                << "  Memory allocated:   " << std::setiosflags(std::ios::fixed)
                                            << std::setprecision(2)
                                            << std::setw(10)
                                            << memoryAllocData << " MB" << std::endl
                << "  Memory used:        " << std::setiosflags(std::ios::fixed)
                                            << std::setprecision(2)
                                            << std::setw(10)
                                            << memoryUsedData << " MB" << std::endl
                << "---------------------------------------------" << std::endl
                << "Adjoint vector"                                << std::endl
                << "---------------------------------------------" << std::endl
                << "  Number of Adjoints: " << std::setw(10) << nAdjoints << std::endl
                << "  Memory allocated:   " << std::setiosflags(std::ios::fixed)
                                            << std::setprecision(2)
                                            << std::setw(10)
                                            << memoryAdjoints << " MB" << std::endl
                << "---------------------------------------------" << std::endl
                << "Indices"                               << std::endl
                << "---------------------------------------------" << std::endl
                << "  Max. live indices:   " << std::setw(10) << maximumGlobalIndex << std::endl
                << "  Cur. live indices:   " << std::setw(10) << currentLiveIndices << std::endl
                << "  Indices stored:      " << std::setw(10) << storedIndices << std::endl
                << "  Memory allocated:    " << std::setiosflags(std::ios::fixed)
                                             << std::setprecision(2)
                                             << std::setw(10)
                                             << memoryAllocatedIndices << " MB" << std::endl
                << "  Memory used:         " << std::setiosflags(std::ios::fixed)
                                             << std::setprecision(2)
                                             << std::setw(10)
                                             << memoryStoredIndices << " MB" << std::endl
                << "-------------------------------------" << std::endl
                << "External functions  "                  << std::endl
                << "-------------------------------------" << std::endl
                << "  Total Number:        " << std::setw(10) << nExternalFunc << std::endl
                << std::endl;

    }

  };
}
