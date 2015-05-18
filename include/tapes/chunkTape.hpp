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
 * Authors: TODO
 */

#pragma once

#include <iostream>
#include <cstddef>
#include <tuple>

#include "activeReal.hpp"
#include "chunk.hpp"
#include "chunkVector.hpp"
#include "reverseTapeInterface.hpp"
#include "externalFunctions.hpp"

namespace codi {

  template <typename IndexType>
  struct ExpressionCounter {

    typedef IndexType Position;

    IndexType count;

    inline Position getPosition() {
      return count;
    }

    inline void reset(const Position& pos) {
      count = pos;
    }
  };

  template <typename Real, typename IndexType>
  struct ChunkTapeTypes {
    typedef Chunk2< Real, IndexType> DataChunk;
    typedef ChunkVector<DataChunk, ExpressionCounter<IndexType> > DataChunkVector;

    typedef Chunk1<OperationInt> OperatorChunk;
    typedef ChunkVector<OperatorChunk, DataChunkVector> OperatorChunkVector;

    typedef Chunk2<ExternalFunction,typename OperatorChunkVector::Position> ExternalFunctionChunk;
    typedef ChunkVector<ExternalFunctionChunk, OperatorChunkVector> ExternalFunctionChunkVector;

    typedef typename ExternalFunctionChunkVector::Position Position;

  };

  template <typename Real, typename IndexType>
  class ChunkTape : public ReverseTapeInterface<Real, IndexType, ChunkTape<Real, IndexType>, typename ChunkTapeTypes<Real, IndexType>::Position > {
  public:

    typedef typename ChunkTapeTypes<Real, IndexType>::DataChunk DataChunk;
    typedef typename ChunkTapeTypes<Real, IndexType>::DataChunkVector DataChunkVector;

    typedef typename ChunkTapeTypes<Real, IndexType>::OperatorChunk OperatorChunk;
    typedef typename ChunkTapeTypes<Real, IndexType>::OperatorChunkVector OperatorChunkVector;

    typedef typename ChunkTapeTypes<Real, IndexType>::ExternalFunctionChunk ExternalFunctionChunk;
    typedef typename ChunkTapeTypes<Real, IndexType>::ExternalFunctionChunkVector ExternalFunctionChunkVector;

    typedef typename ChunkTapeTypes<Real, IndexType>::Position Position;


  private:

    ExpressionCounter<IndexType> expressionCount;
    DataChunkVector data;
    OperatorChunkVector operators;
    ExternalFunctionChunkVector externalFunctions;
    Real* adjoints;
    IndexType adjointsSize;

    bool active;

  public:
    ChunkTape() :
      expressionCount(),
      data(DefaultChunkSize, expressionCount),
      operators(DefaultChunkSize, data),
      externalFunctions(1000, operators),
      adjoints(NULL),
      active(false){
    }

    void setDataChunkSize(const size_t& dataChunkSize) {
      data.setChunkSize(dataChunkSize);
    }

    void setOperatorChunkSize(const size_t& opChunkSize) {
      operators.setChunkSize(opChunkSize);
    }

    void setExternalFunctionChunkSize(const size_t& extChunkSize) {
      externalFunctions.setChunkSize(extChunkSize);
    }

    void resize(const size_t& dataSize, const size_t& opSize) {
      data.resize(dataSize);
      operators.resize(opSize);
    }

    void resizeAdjoints(const IndexType& size) {
      IndexType oldSize = adjointsSize;
      adjointsSize = size;

      adjoints = (Real*)realloc(adjoints, sizeof(Real) * (size_t)adjointsSize);

      for(IndexType i = oldSize; i < adjointsSize; ++i) {
        adjoints[i] = 0.0;
      }
    }

    void allocateAdjoints() {
      resizeAdjoints(expressionCount + 1);
    }

    template<typename Rhs>
    inline void store(Real& lhsValue, IndexType& lhsIndex, const Rhs& rhs) {
      Real gradient; /* This value will not be used */

      ENABLE_CHECK (OptTapeActivity, active){
        data.reserveItems(ExpressionTraits<Rhs>::maxActiveVariables);
        operators.reserveItems(1); // operators needs a reserve bevor the data items for the operator are pushed
        /* first store the size of the current stack position and evaluate the
         rhs expression. If there was an active variable on the rhs, update
         the index of the lhs */
        size_t startSize = data.getChunkPosition();
        rhs.calcGradient(gradient);
        size_t activeVariables = data.getChunkPosition() - startSize;
        if(0 == activeVariables) {
          lhsIndex = 0;
        } else {
          operators.setDataAndMove(std::make_tuple((OperationInt)activeVariables));
          lhsIndex = ++expressionCount.count;
        }
//      } else {
//        lhsIndex = 0;
      }
      /* now set the value of the lhs */
      lhsValue = rhs.getValue();
    }

    inline void store(Real& value, IndexType& lhsIndex, const ActiveReal<Real, ChunkTape<Real, IndexType> >& rhs) {
      ENABLE_CHECK (OptTapeActivity, active){
        lhsIndex = rhs.getGradientData();
//      } else {
//        lhsIndex = 0;
      }
      value = rhs.getValue();
    }

    inline void store(Real& value, IndexType& lhsIndex, const typename TypeTraits<Real>::PassiveReal& rhs) {
      ENABLE_CHECK (OptTapeActivity, active){
        lhsIndex = 0;
      }
      value = rhs;
    }

    inline void pushJacobi(Real& /*gradient*/, const Real& /*value*/, const IndexType& index) {
      if(0 != index) {
        data.setDataAndMove(std::make_tuple(1.0, index));
      }
    }

    inline void pushJacobi(Real& /*gradient*/, const Real& jacobi, const Real& /*value*/, const IndexType& index) {
      if(0 != index) {
        ENABLE_CHECK(OptJacobiIsZero, 0.0 != jacobi) {
          data.setDataAndMove(std::make_tuple(jacobi, index));
        }
      }
    }

    inline void initGradientData(Real& /*value*/, IndexType& index) {
      index = 0;
    }

    inline void destroyGradientData(Real& /*value*/, IndexType& index) {
      /* nothing to do */
    }


    void setGradient(IndexType& index, const Real& gradient) {
      if(0 != index) {
        this->gradient(index) = gradient;
      }
    }

    inline Real getGradient(const IndexType& index) const {
      if(adjointsSize <= index) {
        return Real();
      } else {
        return adjoints[index];
      }
    }

    inline Real& gradient(IndexType& index) {
      if(adjointsSize <= index) {
        resizeAdjoints(index + 1);
      }

      return adjoints[index];
    }

    inline Position getPosition() {
      return externalFunctions.getPosition();
    }

    inline void clearAdjoints(){
      for(IndexType i = 0; i <= expressionCount.count; ++i) {
        adjoints[i] = 0.0;
      }
    }

    inline void clearAdjoints(const Position& start, const Position& end){
      for(IndexType i = start.inner.inner.inner; i <= end.inner.inner.inner; ++i) {
        adjoints[i] = 0.0;
      }
    }

    inline void reset(const Position& pos) {
      for(IndexType i = pos.inner.inner.inner; i <= expressionCount.count; ++i) {
        adjoints[i] = 0.0;
      }

      externalFunctions.forEach(externalFunctions.getPosition(), pos, popExternalFunction);

      // reset will be done iterativly through the vectors
      externalFunctions.reset(pos);
    }

    inline void reset() {
      reset(Position());
    }

    inline void evaluate(const size_t& startAdjPos, const size_t& endAdjPos, size_t& opPos, OperationInt* &operators, size_t& dataPos, Real* &jacobies, IndexType* &indices) {
      size_t adjPos = startAdjPos;

      while(adjPos > endAdjPos) {
        const Real& adj = adjoints[adjPos];
        --adjPos;
        --opPos;
        const OperationInt& activeVariables = operators[opPos];
        ENABLE_CHECK(OptZeroAdjoint, adj != 0){
          for(OperationInt curVar = 0; curVar < activeVariables; ++curVar) {
            --dataPos;
            adjoints[indices[dataPos]] += adj * jacobies[dataPos];

          }
        } else {
          dataPos -= activeVariables;
        }
      }
    }

    inline void evaluateOp(const typename OperatorChunkVector::Position& start, const typename OperatorChunkVector::Position& end) {
      OperationInt* operatorData;
      size_t dataPos = start.data;
      typename DataChunkVector::Position curInnerPos = start.inner;
      for(size_t curChunk = start.chunk; curChunk > end.chunk; --curChunk) {
        std::tie(operatorData) = operators.getDataAtPosition(curChunk, 0);

        typename DataChunkVector::Position endInnerPos = operators.getInnerPosition(curChunk);
        evaluate(curInnerPos, endInnerPos, dataPos, operatorData);

        curInnerPos = endInnerPos;

        dataPos = operators.getChunkUsedData(curChunk - 1);
      }

      // Iterate over the reminder also covers the case if the start chunk and end chunk are the same
      std::tie(operatorData) = operators.getDataAtPosition(end.chunk, 0);
      evaluate(curInnerPos, end.inner, dataPos, operatorData);
    }

    inline void evaluate(const typename DataChunkVector::Position& start, const typename DataChunkVector::Position& end, size_t& opPos, OperationInt* &operatorData) {
      Real* jacobiData;
      IndexType* indexData;
      size_t dataPos = start.data;
      typename ExpressionCounter<IndexType>::Position curInnerPos = start.inner;
      for(size_t curChunk = start.chunk; curChunk > end.chunk; --curChunk) {
        std::tie(jacobiData, indexData) = data.getDataAtPosition(curChunk, 0);

        typename ExpressionCounter<IndexType>::Position endInnerPos = data.getInnerPosition(curChunk);
        evaluate(curInnerPos, endInnerPos, opPos, operatorData, dataPos, jacobiData, indexData);

        curInnerPos = endInnerPos;

        dataPos = data.getChunkUsedData(curChunk - 1);
      }

      // Iterate over the reminder also covers the case if the start chunk and end chunk are the same
      std::tie(jacobiData, indexData) = data.getDataAtPosition(end.chunk, 0);
      evaluate(curInnerPos, end.inner, opPos, operatorData, dataPos, jacobiData, indexData);
    }

    void evaluate(const Position& start, const Position& end) {
      if(adjointsSize <= expressionCount.count) {
        resizeAdjoints(expressionCount.count + 1);
      }

      evaluateExtFunc(start, end);
    }

    struct ExtFuncEvaluator {
      typename OperatorChunkVector::Position curInnerPos;
      ExternalFunction* extFunc;
      typename OperatorChunkVector::Position* endInnerPos;

      ChunkTape<Real, IndexType>& tape;

      ExtFuncEvaluator(typename OperatorChunkVector::Position curInnerPos, ChunkTape<Real, IndexType>& tape) :
        curInnerPos(curInnerPos),
        extFunc(NULL),
        endInnerPos(NULL),
        tape(tape){}

      void operator () (typename ExternalFunctionChunk::DataPointer& data) {
        std::tie(extFunc, endInnerPos) = data;

        // always evaluate the stack to the point of the external function
        tape.evaluateOp(curInnerPos, *endInnerPos);

        if(extFunc->func != NULL){
          extFunc->func(extFunc->checkpoint);
        }

        curInnerPos = *endInnerPos;
      }
    };

    void evaluateExtFunc(const typename ExternalFunctionChunkVector::Position& start, const typename ExternalFunctionChunkVector::Position &end){
      ExtFuncEvaluator evaluator(start.inner, *this);

      externalFunctions.forEach(start, end, evaluator);

      evaluateOp(evaluator.curInnerPos, end.inner);
    }

    void evaluate() {
      evaluate(getPosition(), Position());
    }

    void pushExternalFunctionHandle(ExternalFunction::CallFunction extFunc, void* checkpoint, ExternalFunction::DeleteFunction delCheckpoint){
      externalFunctions.reserveItems(1);
      ExternalFunction function(extFunc, checkpoint, delCheckpoint);
      externalFunctions.setDataAndMove(std::make_tuple(function, operators.getPosition()));
    }

    template<typename Data>
    void pushExternalFunction(typename ExternalFunctionDataHelper<Data>::CallFunction extFunc, Data* checkpoint, typename ExternalFunctionDataHelper<Data>::DeleteFunction delCheckpoint){
      ExternalFunctionDataHelper<Data>* functionHelper = new ExternalFunctionDataHelper<Data>(extFunc, checkpoint, delCheckpoint);
      pushExternalFunctionHandle( ExternalFunctionDataHelper<Data>::callFunction, functionHelper, ExternalFunctionDataHelper<Data>::deleteFunction);
    }

    static void popExternalFunction(typename ExternalFunctionChunk::DataPointer& extFunction) {
      /* we just need to call the delete function */
      std::get<0>(extFunction)->deleteData();
    }

    inline void registerInput(ActiveReal<Real, ChunkTape<Real, IndexType> >& value) {
      operators.reserveItems(1);
      operators.setDataAndMove(std::make_tuple(0));

      value.getGradientData() = ++expressionCount.count;
    }

    inline void registerOutput(ActiveReal<Real, ChunkTape<Real, IndexType> >& /*value*/) {
      /* do nothing */
    }
    inline void setActive(){
      active = true;
    }

    inline void setPassive(){
      active = false;
    }

    inline bool isActive(){
      return active;
    }

  };
}
