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
#include "tapeInterface.hpp"

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
  class ChunkTape : public TapeInterface<Real, IndexType> {
  public:

    typedef Chunk2< Real, IndexType> DataChunk;
    typedef ChunkVector<DataChunk, ExpressionCounter<IndexType> > DataChunkVector;

    typedef Chunk1<IndexType> OperatorChunk;
    typedef ChunkVector<OperatorChunk, DataChunkVector> OperatorChunkVector;

    typedef typename OperatorChunkVector::Position Position;

  private:

    ExpressionCounter<IndexType> expressionCount;
    DataChunkVector data;
    OperatorChunkVector operators;
    Real* adjoints;
    size_t adjointsSize;

    bool active;

  public:
    ChunkTape() :
      expressionCount(),
      data(0, expressionCount),
      operators(0, data),
      adjoints(NULL),
      active(false){
    }

    void resize(const size_t& dataChunkSize, const size_t& opChunkSize) {
      data.setChunkSize(dataChunkSize);
      operators.setChunkSize(opChunkSize);
    }

    void resizeAdjoints(const size_t& size) {
      adjointsSize = size;

      adjoints = (Real*)realloc(adjoints, sizeof(Real) * adjointsSize);
    }

    void allocateAdjoints() {
      resizeAdjoints(expressionCount + 1);
    }

    template<typename Rhs>
    inline void store(Real& lhsValue, IndexType& lhsIndex, const Rhs& rhs) {
      Real gradient; /* This value will not be used */

      if (active){
        data.reserveItems(ExpressionTraits<Rhs>::maxActiveVariables);
        operators.reserveItems(1); // operators needs a reserve bevor the data items for the operator are pushed
        /* first store the size of the current stack position and evaluate the
         rhs expression. If there was an active variable on the rhs, update
         the index of the lhs */
        IndexType startSize = data.getChunkPosition();
        rhs.calcGradient(gradient);
        IndexType activeVariables = data.getChunkPosition() - startSize;
        if(0 == activeVariables) {
          lhsIndex = 0;
        } else {
          operators.setDataAndMove(std::make_tuple(activeVariables));
          lhsIndex = ++expressionCount.count;
        }
      }
      /* now set the value of the lhs */
      lhsValue = rhs.getValue();
    }

    inline void store(Real& value, IndexType& lhsIndex, const ActiveReal<Real, SimpleTape<Real, IndexType> >& rhs) {
      if (active){
        lhsIndex = rhs.getGradientData();
      }
      value = rhs.getValue();
    }

    inline void store(Real& value, IndexType& lhsIndex, const typename TypeTraits<Real>::PassiveReal& rhs) {
      if (active){
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
        data.setDataAndMove(std::make_tuple(jacobi, index));
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
      return operators.getPosition();
    }

    inline void clearAdjoints(){
      for(size_t i = 0; i <= operators.getUsedSize(); ++i) {
        adjoints.data[i] = 0.0;
      }
    }

    inline void reset(const Position& pos) {
      for(size_t i = pos.inner.inner; i <= expressionCount.count; ++i) {
        adjoints[i] = 0.0;
      }

      // reset will be done iterativly through the vectors
      operators.reset(pos);
    }

    inline void reset() {
      reset(Position());
    }

    inline void evaluate(Real* startAdj, Real* endAdj, IndexType* &operators, Real* &jacobies, IndexType* &indices) {
      Real* curAdj = startAdj;

      while(curAdj != endAdj) {
        const Real& adj = *curAdj;
        curAdj--;  // move to next adjoint in array
        operators--; // move to next operator in array
        if (adj != 0){
          for(IndexType curVar = 0; curVar < *operators; ++curVar) {
            indices--;  // move to next index in array
            jacobies--; // move to next jacobi in array
            adjoints[*indices] += adj * *jacobies;

          }
        }else {
          indices -= *operators;
          jacobies -= *operators;
        }
      }
    }

    inline void evaluateOp(const typename OperatorChunkVector::Position& start, const typename OperatorChunkVector::Position& end) {
      IndexType* operatorPos;
      size_t dataPos = start.data;
      typename DataChunkVector::Position curInnerPos = start.inner;
      for(size_t curChunk = start.chunk; curChunk > end.chunk; --curChunk) {
        std::tie(operatorPos) = operators.getDataAtPosition(curChunk, dataPos);

        typename DataChunkVector::Position endInnerPos = operators.getInnerPosition(curChunk);
        evaluate(curInnerPos, endInnerPos, operatorPos);

        curInnerPos = endInnerPos;

        dataPos = operators.getChunkUsedData(curChunk - 1);
      }

      // Iterate over the reminder also covers the case if the start chunk and end chunk are the same
      std::tie(operatorPos) = operators.getDataAtPosition(end.chunk, dataPos);
      evaluate(curInnerPos, end.inner, operatorPos);
    }

    inline void evaluate(const typename DataChunkVector::Position& start, const typename DataChunkVector::Position& end, IndexType* &operatorPos) {
      Real* jacobiPos;
      IndexType* indexPos;
      size_t dataPos = start.data;
      typename ExpressionCounter<IndexType>::Position curInnerPos = start.inner;
      for(size_t curChunk = start.chunk; curChunk > end.chunk; --curChunk) {
        std::tie(jacobiPos, indexPos) = data.getDataAtPosition(curChunk, dataPos);

        typename ExpressionCounter<IndexType>::Position endInnerPos = data.getInnerPosition(curChunk);
        evaluate(&adjoints[curInnerPos], &adjoints[endInnerPos], operatorPos, jacobiPos, indexPos);

        curInnerPos = endInnerPos;

        dataPos = data.getChunkUsedData(curChunk - 1);
      }

      // Iterate over the reminder also covers the case if the start chunk and end chunk are the same
      std::tie(jacobiPos, indexPos) = data.getDataAtPosition(end.chunk, dataPos);
      evaluate(&adjoints[curInnerPos], &adjoints[end.inner], operatorPos, jacobiPos, indexPos);
    }

    void evaluate(const Position& start, const Position& end) {
      if(adjointsSize <= expressionCount.count) {
        resizeAdjoints(expressionCount.count + 1);
      }

      evaluateOp(start, end);
    }

    void evaluate() {
      evaluate(getPosition(), Position());
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
