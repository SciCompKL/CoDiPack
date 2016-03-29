/**
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
#include <tuple>

#include "../activeReal.hpp"
#include "../expressionHandle.hpp"
#include "chunk.hpp"
#include "indices/linearIndexHandler.hpp"
#include "reverseTapeInterface.hpp"

namespace codi {

  /**
   * @brief Position for the simple tape.
   */
  struct SimplePrimalValueTapePosition {
    /** @brief The current statement recorded on the tape. */
    size_t stmt;
    /** @brief The current jacobi data recorded on the tape. */
    size_t data;
    /** @brief The current passive value data recorded on the tape. */
    size_t passiveData;
    /** @brief The current external function recorded on the tape. */
    size_t extFunc;

    /**
     * @brief Simple constructor for convenience.
     * @param[in]    stmt  The current statement recorded on the tape.
     * @param[in]    data  The current jacobi recorded on the tape.
     * @param[in] extFunc  The current external function recorded on the tape.
     */
    SimplePrimalValueTapePosition(const size_t& stmt, const size_t& data, const size_t& passiveData, const size_t& extFunc) :
      stmt(stmt),
      data(data),
      passiveData(passiveData),
      extFunc(extFunc) {}

    SimplePrimalValueTapePosition() :
      stmt(0),
      data(0),
      passiveData(0),
      extFunc(0) {}
  };

  template<typename IndexType, size_t n>
  struct PassiveDataHelper {
    size_t pos;
    IndexType indices[n];

    PassiveDataHelper() : pos(0) {}

    inline void push(const IndexType& index) {
      indices[pos++] = index;
    }

    inline void reset() {
      pos = 0;
    }
  };

  /**
   * @brief A tape with a simple implementation and no bounds checking.
   *
   * The SimplePrimalValueTape implements a fully featured ReverseTapeInterface in a
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
  template <typename RealType, typename IndexType>
  class SimplePrimalValueTape : public ReverseTapeInterface<RealType, IndexType, RealType, SimplePrimalValueTape<RealType, IndexType>, SimplePrimalValueTapePosition > {
  public:

    /** @brief The type for the primal values. */
    typedef RealType Real;
    /** @brief The type for the adjoint values. */
    typedef Real GradientValue;
    /** @brief The type for the index handler. */
    typedef LinearIndexHandler<IndexType> IndexHandler;

    /** @brief The type for the indices that are used for the identification of the adjoint variables. */
    //typedef typename IndexHandler::IndexType IndexType;
    /** @brief The gradient data is just the index type. */
    typedef IndexType GradientData;

    typedef typename TypeTraits<Real>::PassiveReal PassiveReal;

    #define TAPE_NAME SimplePrimalValueTape

    #define POSITION_TYPE SimplePrimalValueTapePosition
    #define INDEX_HANDLER_TYPE IndexHandler
    #define RESET_FUNCTION_NAME resetInt
    #define EVALUATE_FUNCTION_NAME evaluateInt
    #include "modules/tapeBaseModule.tpp"

    template<typename AdjointData>
    static void inputHandleFunc(const Real& seed, const IndexType* indices, const PassiveReal* passiveValues, const Real* primalValues, Real* adjointValues) {}
    const static ExpressionHandle<Real*, Real, IndexType> InputHandle;

    Chunk1<IndexType> data;
    Chunk1<PassiveReal> passiveData;
    Chunk1<const ExpressionHandle<Real*, Real, IndexType>* > statements;
    /**
     * @brief The external function data and the position where the external function has been inserted.
     */
    Chunk2<ExternalFunction, Position> externalFunctions;

    Real* primals;
    IndexType primalsSize;
    IndexType primalsIncr;

    PassiveDataHelper<IndexType, 256> passiveDataHelper;

  public:
    /**
     * @brief Creates a tape with the size of zero for the data, statements and external functions.
     */
    SimplePrimalValueTape() :
      /* defined in tapeBaseModule */indexHandler(),
      /* defined in tapeBaseModule */adjoints(NULL),
      /* defined in tapeBaseModule */adjointsSize(0),
      /* defined in tapeBaseModule */active(false),
      data(0),
      passiveData(0),
      statements(0),
      externalFunctions(0),
      primals(NULL),
      primalsSize(0),
      primalsIncr(DefaultSmallChunkSize) {}

    /**
     * @brief Helper function: Sets the primal vector to a new size.
     *
     * @param[in] size The new size for the primal vector.
     */
    void resizePrimals(const IndexType& size) {
      IndexType oldSize = primalsSize;
      primalsSize = size;

      primals = (Real*)realloc(adjoints, sizeof(Real) * (size_t)primalsSize);
    }

    inline void checkPrimalsSize() {
      //if(primalsSize <= indexHandler.getMaximumGlobalIndex()) {
      if(primalsSize <= statements.getUsedSize()) {
        IndexType newSize = 1 + (statements.getUsedSize() + 1) / primalsIncr;
        newSize = newSize * primalsIncr;
        resizePrimals(newSize);
      }
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
      for(IndexType i = end.stmt; i <= start.stmt; ++i) {
        adjoints[i] = GradientValue();
      }
    }

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

    void setPassiveDataSize(const size_t& passiveDataSize) {
      passiveData.resize(passiveDataSize);
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
     * @brief Return the number of passive data entries.
     * @return The number of passive data entries.
     */
    size_t getUsedPassiveDataSize() {
      return passiveData.getUsedSize();
    }

    /**
     * @brief Set the size of the jacobi and statement data and the adjoint vector.
     * @param[in] dataSize  The new size of the jacobi vector.
     * @param[in] stmtSize  The new size of the statement vector.
     */
    void resize(const size_t& dataSize, const size_t& stmtSize) {
      data.resize(dataSize);
      statements.resize(stmtSize);

      resizePrimals(stmtSize + 1);
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

      ENABLE_CHECK(OptTapeActivity, active){
        codiAssert(ExpressionTraits<Rhs>::maxActiveVariables <= data.getUnusedSize());
        codiAssert(ExpressionTraits<Rhs>::maxPassiveVariables <= passiveData.getUnusedSize());
        size_t dataSize = data.getUsedSize();
        size_t passiveDataSize = passiveData.getUsedSize();
        CODI_UNUSED(dataSize);  /* needed to avoid unused variable when the assersts are not enabled. */
        CODI_UNUSED(passiveDataSize);  /* needed to avoid unused variable when the assersts are not enabled. */
        rhs.pushIndices(this);
        rhs.pushPassive(this);
        codiAssert(ExpressionTraits<Rhs>::maxActiveVariables == data.getUsedSize() - dataSize);
        codiAssert(ExpressionTraits<Rhs>::maxPassiveVariables == passiveData.getUsedSize() - passiveDataSize);
        codiAssert(statements.getUsedSize() < statements.size);
        statements.setDataAndMove(ExpressionHandleStore<Real*, Real, IndexType, Rhs>::getHandle());
        checkPrimalsSize();
        primals[statements.getUsedSize()] = rhs.getValue();
        lhsIndex = statements.getUsedSize();

        // clear the generated temproal indices
        if(0 != passiveDataHelper.pos) {
          for(size_t i = 0; i < passiveDataHelper.pos; ++i) {
            destroyGradientData(primals[passiveDataHelper.indices[i]], passiveDataHelper.indices[i]);
          }
          passiveDataHelper.reset();
        }
      } else {
        lhsIndex = 0;
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
    inline void store(Real& lhsValue, IndexType& lhsIndex, const ActiveReal<SimplePrimalValueTape<Real, IndexType> >& rhs) {
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
      ENABLE_CHECK(OptTapeActivity, active){
        codiAssert(statements.getUsedSize() < statements.size);

        statements.setDataAndMove(&InputHandle);
        lhsIndex = statements.getUsedSize();

        checkPrimalsSize();
        primals[statements.getUsedSize()] = rhs;
      } else {
        lhsIndex = 0;
      }

      lhsValue = rhs;
    }

    inline void pushPassive(const PassiveReal& value) {
      codiAssert(passiveData.getUsedSize() < passiveData.size);

      passiveData.setDataAndMove(value);
    }

    inline void pushIndices(const Real& value, const IndexType& index) {
      codiAssert(data.getUsedSize() < data.size);

      if(0 == index) {
        codiAssert(statements.getUsedSize() < statements.size);
        // create temporary index
        statements.setDataAndMove(&InputHandle);
        IndexType tempIndex = statements.getUsedSize();
        checkPrimalsSize();
        primals[tempIndex] = value;
        data.setDataAndMove(tempIndex);

        passiveDataHelper.push(index);
      } else {
        data.setDataAndMove(index);
      }
    }

    /**
     * @brief Stores the jacobi with the value 1.0 on the tape if the index is active.
     *
     * @param[in] gradient Not used in this implementation.
     * @param[in]    value Not used in this implementation.
     * @param[in]    index Used to check if the variable is active.
     */
    template<typename Data>
    inline void pushJacobi(Data& passiveDataHelper, const Real& value, const IndexType& index) {
      CODI_UNUSED(passiveDataHelper);
      CODI_UNUSED(value);
      CODI_UNUSED(index);

      codiAssert(false || "Should not be called.");
    }

    /**
     * @brief Stores the jacobi on the tape if the index is active.
     *
     * @param[in] gradient Not used in this implementation.
     * @param[in]   jacobi Stored on the tape if the variable is active.
     * @param[in]    value Not used in this implementation.
     * @param[in]    index Used to check if the variable is active.
     */
    template<typename Data>
    inline void pushJacobi(Data& passiveDataHelper, const Real& jacobi, const Real& value, const IndexType& index) {
      CODI_UNUSED(passiveDataHelper);
      CODI_UNUSED(jacobi);
      CODI_UNUSED(value);
      CODI_UNUSED(index);

      codiAssert(false || "Should not be called.");
    }

    /**
     * @brief Get the current position of the tape.
     *
     * The position can be used to reset the tape to that position or to
     * evaluate only parts of the tape.
     * @return The current position of the tape.
     */
    inline Position getPosition() const {
      return Position(statements.getUsedSize(), data.getUsedSize(), passiveData.getUsedSize(), externalFunctions.getUsedSize());
    }

    /**
     * @brief Reset the tape to the given position.
     *
     * @param[in] pos Reset the state of the tape to the given position.
     */
    inline void resetInt(const Position& pos) {
      codiAssert(pos.stmt < statements.size);
      codiAssert(pos.data < data.size);
      codiAssert(pos.passiveData < passiveData.size);
      codiAssert(pos.extFunc < externalFunctions.size);

      for(size_t i = pos.stmt; i <= primalsSize; ++i) {
        primals[i] = 0.0;
      }

      for(size_t i = pos.extFunc; i < externalFunctions.getUsedSize(); ++i) {
        externalFunctions.data1[i].deleteData();
      }

      statements.setUsedSize(pos.stmt);
      data.setUsedSize(pos.data);
      passiveData.setUsedSize(pos.passiveData);
      externalFunctions.setUsedSize(pos.extFunc);
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
        const Real& adj = adjoints[curPos.stmt];
        --curPos.stmt;
        const ExpressionHandle<Real*, Real, IndexType>* exprHandle = statements.data[curPos.stmt];
        curPos.data -= exprHandle->maxActiveVariables;
        curPos.passiveData -= exprHandle->maxPassiveVariables;
        ENABLE_CHECK(OptZeroAdjoint, adj != 0.0){

          const IndexType* indices = &data.data[curPos.data];
          const PassiveReal* passiveValues = &passiveData.data[curPos.passiveData];
          exprHandle->adjointFunc(adj, indices, passiveValues, primals, adjoints);
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
    inline void evaluateInt(const Position& start, const Position& end) {
      codiAssert(start.data >= end.data);
      codiAssert(start.stmt >= end.stmt);
      codiAssert(start.extFunc >= end.extFunc);

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
     * @brief Register a variable as an active variable.
     *
     * The index of the variable is set to a non zero number.
     * @param[inout] value The value which will be marked as an active variable.
     */
    inline void registerInput(ActiveReal<SimplePrimalValueTape<Real, IndexType> >& value) {
      codiAssert(statements.getUsedSize() < statements.size);

      statements.setDataAndMove(&InputHandle);
      value.getGradientData() = statements.getUsedSize();

      checkPrimalsSize();
      primals[statements.getUsedSize()] = value.getValue();
    }

    /**
     * @brief Not needed in this implementation.
     *
     * @param[in] value Not used.
     */
    inline void registerOutput(ActiveReal<SimplePrimalValueTape<Real, IndexType> >& value) {
      value = 1.0 * value;
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

    template<typename Stream = std::ostream>
    void printStatistics(Stream& out = std::cout) const {

      const std::string hLine = "-------------------------------------\n";

      size_t nChunksData  = 1;
      size_t totalData    =  data.getUsedSize();
      double  memoryUsedData = (double)totalData*(double)(sizeof(IndexType))* BYTE_TO_MB;
      double  memoryAllocData= (double)data.getSize() *(double)(sizeof(IndexType))* BYTE_TO_MB;

      size_t nChunksStmts  = 1;
      size_t totalStmts    = statements.getUsedSize();

      double  memoryUsedStmts = (double)totalStmts*(double)sizeof(const ExpressionHandle<Real*, Real, IndexType>*)* BYTE_TO_MB;
      double  memoryAllocStmts= (double)statements.getSize()*(double)sizeof(const ExpressionHandle<Real*, Real, IndexType>*)* BYTE_TO_MB;

      size_t nChunksPassive  = 1;
      size_t totalPassive    =  passiveData.getUsedSize();
      double  memoryUsedPassive = (double)totalPassive*(double)(sizeof(PassiveReal))* BYTE_TO_MB;
      double  memoryAllocPassive = (double)passiveData.getSize() *(double)(sizeof(PassiveReal))* BYTE_TO_MB;

      out << hLine
          << "CoDi Tape Statistics ( simple primal value tape)\n";
      printTapeBaseStatistics(out, hLine);
      out << hLine
          << "Statements \n"
          << hLine
          << "  Number of Chunks: " << std::setw(10) << nChunksStmts << "\n"
          << "  Total Number:     " << std::setw(10) << totalStmts   << "\n"
          << "  Memory allocated: " << std::setiosflags(std::ios::fixed)
                                    << std::setprecision(2)
                                    << std::setw(10)
                                    << memoryAllocStmts << " MB" << "\n"
          << "  Memory used:      " << std::setiosflags(std::ios::fixed)
                                    << std::setprecision(2)
                                    << std::setw(10)
                                    << memoryUsedStmts << " MB" << "\n";
      out << hLine
          << "Data entries \n"
          << hLine
          << "  Number of Chunks: " << std::setw(10) << nChunksData << "\n"
          << "  Total Number:     " << std::setw(10) << totalData   << "\n"
          << "  Memory allocated: " << std::setiosflags(std::ios::fixed)
                                    << std::setprecision(2)
                                    << std::setw(10)
                                    << memoryAllocData << " MB" << "\n"
          << "  Memory used:      " << std::setiosflags(std::ios::fixed)
                                    << std::setprecision(2)
                                    << std::setw(10)
                                    << memoryUsedData << " MB" << "\n";
      out << hLine
          << "Passive data entries \n"
          << hLine
          << "  Number of Chunks: " << std::setw(10) << nChunksPassive << "\n"
          << "  Total Number:     " << std::setw(10) << totalPassive << "\n"
          << "  Memory allocated: " << std::setiosflags(std::ios::fixed)
                                    << std::setprecision(2)
                                    << std::setw(10)
                                    << memoryAllocPassive << " MB" << "\n"
          << "  Memory used:      " << std::setiosflags(std::ios::fixed)
                                    << std::setprecision(2)
                                    << std::setw(10)
                                    << memoryUsedPassive << " MB" << "\n";


    }

  private:
    /**
     * @brief Private common method to add to the external function stack.
     *
     * @param[in] function The external function structure to push.
     */
    void pushExternalFunctionHandle(const ExternalFunction& function){
      codiAssert(0 != externalFunctions.getUnusedSize());
      externalFunctions.setDataAndMove(function, getPosition());
    }

  };

  template <typename Real, typename IndexType>
  const ExpressionHandle<Real*, Real, IndexType> SimplePrimalValueTape<Real, IndexType>::InputHandle(&SimplePrimalValueTape<Real, IndexType>::inputHandleFunc<Real*>, 0, 0);
}
