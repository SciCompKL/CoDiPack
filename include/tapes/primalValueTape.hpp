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
   * @brief Vector defintion for the SimplePrimalValueTape.
   *
   * The structure defines all vectors as single chunk vectors.
   *
   * See PrimalValueTape for details.
   *
   * @tparam Real  The type for the primal values.
   * @tparam IndexHandler  The index handler for the managing of the indices. It has to be a index handler that assumes index reuse.
   * @tparam GradientValue  The type for the adjoint values. (Default: Same as the primal value.)
   */
  template <typename Real, typename IndexHandler, typename GradientValue = Real>
  struct SimplePrimaValuelTapeTypes {
    /** @brief The type for the primal values. */
    typedef Real RealType;
    /** @brief The handler for the indices. */
    typedef IndexHandler IndexHandlerType;
    /** @brief The type for the adjoint values. */
    typedef GradientValue GradientValueType;

    typedef typename IndexHandler::IndexType IndexType;

    typedef const ExpressionHandle<Real*, Real, IndexType>* HandleType;

    /** @brief The data for each statement. */
    typedef Chunk1<HandleType> StatementChunk;
    /** @brief The chunk vector for the statement data. */
    typedef SingleChunkVector<StatementChunk, IndexHandler> StatementVector;

    /** @brief The data for the indices of each statement */
    typedef Chunk1< typename IndexHandler::IndexType> IndexChunk;
    /** @brief The chunk vector for the index data. */
    typedef SingleChunkVector<IndexChunk, StatementVector> IndexVector;

    /** @brief The data for the passive values of each statement */
    typedef Chunk1< typename TypeTraits<Real>::PassiveReal> PassiveChunk;
    /** @brief The chunk vector for the passive data. */
    typedef SingleChunkVector<PassiveChunk, IndexVector> PassiveVector;

    /** @brief The data for the external functions. */
    typedef Chunk2<ExternalFunction,typename PassiveVector::Position> ExternalFunctionChunk;
    /** @brief The chunk vector for the external  function data. */
    typedef SingleChunkVector<ExternalFunctionChunk, PassiveVector> ExternalFunctionVector;

    /** @brief The position for all the different data vectors. */
    typedef typename ExternalFunctionVector::Position Position;

    constexpr static const char* tapeName = "SimplePrimalValueTape";

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
  template <typename TapeTypes>
  class PrimalValueTape : public ReverseTapeInterface<typename TapeTypes::RealType, typename TapeTypes::IndexHandlerType::IndexType, typename TapeTypes::GradientValueType, PrimalValueTape<TapeTypes>, typename TapeTypes::Position > {
  public:

    /** @brief The type for the primal values. */
    typedef typename TapeTypes::RealType Real;
    /** @brief The type for the adjoint values. */
    typedef typename TapeTypes::GradientValueType GradientValue;
    /** @brief The type for the index handler. */
    typedef typename TapeTypes::IndexHandlerType IndexHandler;

    /** @brief The type for the indices that are used for the identification of the adjoint variables. */
    typedef typename TapeTypes::IndexType IndexType;
    /** @brief The gradient data is just the index type. */
    typedef IndexType GradientData;

    typedef typename TypeTraits<Real>::PassiveReal PassiveReal;
    typedef typename TapeTypes::HandleType Handle;

    typedef typename TapeTypes::PassiveVector PassiveVector;
    typedef typename TapeTypes::IndexVector IndexVector;
    typedef typename TapeTypes::StatementVector StatementVector;

    typedef typename TapeTypes::PassiveVector::Position PassivePosition;
    typedef typename TapeTypes::IndexVector::Position IndexPosition;
    typedef typename TapeTypes::StatementVector::Position StmtPosition;

    #define TAPE_NAME PrimalValueTape

    #define POSITION_TYPE typename TapeTypes::Position
    #define INDEX_HANDLER_TYPE IndexHandler
    #define RESET_FUNCTION_NAME resetInt
    #define EVALUATE_FUNCTION_NAME evaluateExtFunc
    #include "modules/tapeBaseModule.tpp"

    typename TapeTypes::StatementVector stmtVector;
    typename TapeTypes::IndexVector indexVector;
    typename TapeTypes::PassiveVector passiveVector;

    #define CHILD_VECTOR_TYPE PassiveVector
    #define CHILD_VECTOR_NAME passiveVector
    #define VECTOR_TYPE typename TapeTypes::ExternalFunctionVector
    #include "modules/externalFunctionsModule.tpp"

    #undef TAPE_NAME

    static void inputHandleFunc(const Real& seed, const IndexType* indices, const PassiveReal* passiveValues, const Real* primalValues, Real* adjointValues) {}
    const static ExpressionHandle<Real*, Real, IndexType> InputHandle;

    Real* primals;
    IndexType primalsSize;
    IndexType primalsIncr;

    PassiveDataHelper<IndexType, 256> passiveDataHelper;

  public:
    /**
     * @brief Creates a tape with the size of zero for the data, statements and external functions.
     */
    PrimalValueTape() :
      /* defined in tapeBaseModule */indexHandler(),
      /* defined in tapeBaseModule */adjoints(NULL),
      /* defined in tapeBaseModule */adjointsSize(0),
      /* defined in tapeBaseModule */active(false),
      stmtVector(DefaultChunkSize, indexHandler),
      indexVector(DefaultChunkSize, stmtVector),
      passiveVector(DefaultChunkSize, indexVector),
      /* defined in externalFunctionsModule */extFuncVector(1000, passiveVector),
      primals(NULL),
      primalsSize(0),
      primalsIncr(DefaultSmallChunkSize){}

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
      if(primalsSize <= indexHandler.getMaximumGlobalIndex()) {
        IndexType newSize = 1 + (indexHandler.getMaximumGlobalIndex() + 1) / primalsIncr;
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
      for(IndexType i = end.inner.inner.inner.inner; i <= start.inner.inner.inner.inner; ++i) {
        adjoints[i] = GradientValue();
      }
    }

    void setPassiveDataSize(const size_t& passiveDataSize) {
      passiveVector.resize(passiveDataSize);
    }

    /**
     * @brief Set the size of the jacobi and statement data and the adjoint vector.
     * @param[in] dataSize  The new size of the jacobi vector.
     * @param[in] stmtSize  The new size of the statement vector.
     */
    void resize(const size_t& dataSize, const size_t& stmtSize) {
      indexVector.resize(dataSize);
      stmtVector.resize(stmtSize);

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
        indexVector.reserveItems(ExpressionTraits<Rhs>::maxActiveVariables);
        size_t indexSize = indexVector.getChunkPosition();
        CODI_UNUSED(indexSize);  /* needed to avoid unused variable when the assersts are not enabled. */
        rhs.pushIndices(this);
        codiAssert(ExpressionTraits<Rhs>::maxActiveVariables == indexVector.getChunkPosition() - indexSize);

        passiveVector.reserveItems(ExpressionTraits<Rhs>::maxPassiveVariables);
        size_t passiveSize = passiveVector.getChunkPosition();
        CODI_UNUSED(passiveSize);  /* needed to avoid unused variable when the assersts are not enabled. */
        rhs.pushPassive(this);
        codiAssert(ExpressionTraits<Rhs>::maxPassiveVariables == passiveVector.getChunkPosition() - passiveSize);

        stmtVector.reserveItems(1);
        stmtVector.setDataAndMove(ExpressionHandleStore<Real*, Real, IndexType, Rhs>::getHandle());
        indexHandler.assignIndex(lhsIndex);

        checkPrimalsSize();
        primals[lhsIndex] = rhs.getValue();

        // clear the generated temporal indices
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
    inline void store(Real& lhsValue, IndexType& lhsIndex, const ActiveReal<PrimalValueTape<TapeTypes> >& rhs) {
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
        // the default behaviour is to activate passive assignments in order to have less passive values
        // in the store for expression routines
        stmtVector.reserveItems(1);
        stmtVector.setDataAndMove(&InputHandle);
        indexHandler.assignIndex(lhsIndex);

        checkPrimalsSize();
        primals[lhsIndex] = rhs;
      } else {
        indexHandler.freeIndex(lhsIndex);
      }

      lhsValue = rhs;
    }

    inline void pushPassive(const PassiveReal& value) {
      passiveVector.setDataAndMove(value);
    }

    inline void pushIndices(const Real& value, const IndexType& index) {
      if(0 == index) {
        // create temporary index
        stmtVector.reserveItems(1);
        stmtVector.setDataAndMove(&InputHandle);

        IndexType tempIndex = 0;
        indexHandler.assignIndex(tempIndex);

        checkPrimalsSize();
        primals[tempIndex] = value;

        passiveDataHelper.push(tempIndex);

        indexVector.setDataAndMove(tempIndex);
      } else {
        indexVector.setDataAndMove(index);
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
      return getExtFuncPosition();
    }

    /**
     * @brief Reset the tape to the given position.
     *
     * @param[in] pos Reset the state of the tape to the given position.
     */
    inline void resetInt(const Position& pos) {
      for(IndexType i = pos.inner.inner.inner.inner + 1; i < primalsSize; ++i) {
        primals[i] = 0.0;
      }

      resetExtFunc(pos);

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
    inline void evaluateStack(const size_t& startAdjPos, const size_t& endAdjPos, size_t& stmtPos, Handle* &statements, size_t& indexPos, IndexType* &indices, size_t& passivePos, PassiveReal* &passives) {
      size_t adjPos = startAdjPos;

      while(adjPos > endAdjPos) {
        const Real& adj = adjoints[adjPos];
        --adjPos;
        --stmtPos;
        Handle exprHandle = statements[stmtPos];
        indexPos -= exprHandle->maxActiveVariables;
        passivePos -= exprHandle->maxPassiveVariables;
        ENABLE_CHECK(OptZeroAdjoint, adj != 0.0){

          exprHandle->adjointFunc(adj, &indices[indexPos], &passives[passivePos], primals, adjoints);
        }
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
    inline void evalStmt(const StmtPosition& start, const StmtPosition& end, Args&&... args) {
      Handle* data;
      size_t dataPos = start.data;
      auto curInnerPos = start.inner;
      for(size_t curChunk = start.chunk; curChunk > end.chunk; --curChunk) {
        stmtVector.getDataAtPosition(curChunk, 0, data);

        auto endInnerPos = stmtVector.getInnerPosition(curChunk);
        evaluateStack(curInnerPos, endInnerPos, dataPos, data, std::forward<Args>(args)...);

        curInnerPos = endInnerPos;

        dataPos = stmtVector.getChunkUsedData(curChunk - 1);
      }

      // Iterate over the reminder also covers the case if the start chunk and end chunk are the same
      stmtVector.getDataAtPosition(end.chunk, 0, data);
      evaluateStack(curInnerPos, end.inner, dataPos, data, std::forward<Args>(args)...);
    }

    template<typename ... Args>
    inline void evalIndices(const IndexPosition& start, const IndexPosition& end, Args&&... args) {
      IndexType* data;
      size_t dataPos = start.data;
      auto curInnerPos = start.inner;
      for(size_t curChunk = start.chunk; curChunk > end.chunk; --curChunk) {
        indexVector.getDataAtPosition(curChunk, 0, data);

        auto endInnerPos = indexVector.getInnerPosition(curChunk);
        evalStmt(curInnerPos, endInnerPos, dataPos, data, std::forward<Args>(args)...);

        curInnerPos = endInnerPos;

        dataPos = indexVector.getChunkUsedData(curChunk - 1);
      }

      // Iterate over the reminder also covers the case if the start chunk and end chunk are the same
      indexVector.getDataAtPosition(end.chunk, 0, data);
      evalStmt(curInnerPos, end.inner, dataPos, data, std::forward<Args>(args)...);
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
    inline void evalExtFuncCallback(const PassivePosition& start, const PassivePosition& end) {
      PassiveReal* data;
      size_t dataPos = start.data;
      auto curInnerPos = start.inner;
      for(size_t curChunk = start.chunk; curChunk > end.chunk; --curChunk) {
        passiveVector.getDataAtPosition(curChunk, 0, data);

        auto endInnerPos = passiveVector.getInnerPosition(curChunk);
        evalIndices(curInnerPos, endInnerPos, dataPos, data);

        curInnerPos = endInnerPos;

        dataPos = passiveVector.getChunkUsedData(curChunk - 1);
      }

      // Iterate over the reminder also covers the case if the start chunk and end chunk are the same
      passiveVector.getDataAtPosition(end.chunk, 0, data);
      evalIndices(curInnerPos, end.inner, dataPos, data);
    }

  public:

    /**
     * @brief Register a variable as an active variable.
     *
     * The index of the variable is set to a non zero number.
     * @param[inout] value The value which will be marked as an active variable.
     */
    inline void registerInput(ActiveReal<PrimalValueTape<TapeTypes> >& value) {
      stmtVector.reserveItems(1);
      stmtVector.setDataAndMove(&InputHandle);
      indexHandler.assignIndex(value.getGradientData());

      checkPrimalsSize();
      primals[value.getGradientData()] = value.getValue();
    }

    /**
     * @brief Not needed in this implementation.
     *
     * @param[in] value Not used.
     */
    inline void registerOutput(ActiveReal<PrimalValueTape<TapeTypes> >& value) {
      value = 1.0 * value;
    }

    template<typename Stream = std::ostream>
    void printStatistics(Stream& out = std::cout) const {

      const std::string hLine = "-------------------------------------\n";

//      size_t nChunksData  = 1;
//      size_t totalData    =  data.getUsedSize();
//      double  memoryUsedData = (double)totalData*(double)(sizeof(IndexType))* BYTE_TO_MB;
//      double  memoryAllocData= (double)data.getSize() *(double)(sizeof(IndexType))* BYTE_TO_MB;
//
//      size_t nChunksStmts  = 1;
//      size_t totalStmts    = statements.getUsedSize();
//
//      double  memoryUsedStmts = (double)totalStmts*(double)sizeof(const ExpressionHandle<Real*, Real, IndexType>*)* BYTE_TO_MB;
//      double  memoryAllocStmts= (double)statements.getSize()*(double)sizeof(const ExpressionHandle<Real*, Real, IndexType>*)* BYTE_TO_MB;
//
//      size_t nChunksPassive  = 1;
//      size_t totalPassive    =  passiveData.getUsedSize();
//      double  memoryUsedPassive = (double)totalPassive*(double)(sizeof(PassiveReal))* BYTE_TO_MB;
//      double  memoryAllocPassive = (double)passiveData.getSize() *(double)(sizeof(PassiveReal))* BYTE_TO_MB;

      out << hLine
          << "CoDi Tape Statistics (simple primal value tape)\n";
      printTapeBaseStatistics(out, hLine);
//      out << hLine
//          << "Statements \n"
//          << hLine
//          << "  Number of Chunks: " << std::setw(10) << nChunksStmts << "\n"
//          << "  Total Number:     " << std::setw(10) << totalStmts   << "\n"
//          << "  Memory allocated: " << std::setiosflags(std::ios::fixed)
//                                    << std::setprecision(2)
//                                    << std::setw(10)
//                                    << memoryAllocStmts << " MB" << "\n"
//          << "  Memory used:      " << std::setiosflags(std::ios::fixed)
//                                    << std::setprecision(2)
//                                    << std::setw(10)
//                                    << memoryUsedStmts << " MB" << "\n";
//      out << hLine
//          << "Data entries \n"
//          << hLine
//          << "  Number of Chunks: " << std::setw(10) << nChunksData << "\n"
//          << "  Total Number:     " << std::setw(10) << totalData   << "\n"
//          << "  Memory allocated: " << std::setiosflags(std::ios::fixed)
//                                    << std::setprecision(2)
//                                    << std::setw(10)
//                                    << memoryAllocData << " MB" << "\n"
//          << "  Memory used:      " << std::setiosflags(std::ios::fixed)
//                                    << std::setprecision(2)
//                                    << std::setw(10)
//                                    << memoryUsedData << " MB" << "\n";
//      out << hLine
//          << "Passive data entries \n"
//          << hLine
//          << "  Number of Chunks: " << std::setw(10) << nChunksPassive << "\n"
//          << "  Total Number:     " << std::setw(10) << totalPassive << "\n"
//          << "  Memory allocated: " << std::setiosflags(std::ios::fixed)
//                                    << std::setprecision(2)
//                                    << std::setw(10)
//                                    << memoryAllocPassive << " MB" << "\n"
//          << "  Memory used:      " << std::setiosflags(std::ios::fixed)
//                                    << std::setprecision(2)
//                                    << std::setw(10)
//                                    << memoryUsedPassive << " MB" << "\n";
      printExternalFunctionStatistics(out, hLine);


    }
  };

  template <typename TapeTypes>
  const ExpressionHandle<typename TapeTypes::RealType*, typename TapeTypes::RealType, typename TapeTypes::IndexType> PrimalValueTape<TapeTypes>::InputHandle(&PrimalValueTape<TapeTypes>::inputHandleFunc, 0, 0);
}
