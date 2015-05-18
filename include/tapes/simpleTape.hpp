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

#include <cstddef>
#include <tuple>

#include "activeReal.hpp"
#include "chunk.hpp"
#include "reverseTapeInterface.hpp"

namespace codi {

  struct SimpleTapePosition {
    size_t stmt;
    size_t data;
    size_t extFunc;

    SimpleTapePosition(const size_t& stmt, const size_t& data, const size_t& extFunc) :
      stmt(stmt),
      data(data),
      extFunc(extFunc) {}
  };

  template <typename Real, typename IndexType>
  class SimpleTape : public ReverseTapeInterface<Real, IndexType, SimpleTape<Real, IndexType>, SimpleTapePosition > {
  public:

    typedef SimpleTapePosition Position;
  private:
    Chunk2<Real, IndexType> data;
    Chunk1<StatementInt> statements;
    Chunk2<ExternalFunction, Position> externalFunctions;
    Chunk1<Real> adjoints;

    bool active;

  public:
    SimpleTape() :
      data(0),
      statements(0),
      externalFunctions(0),
      adjoints(1),
      active(false){}

    void setExternalFunctionChunkSize(const size_t& extChunkSize) {
      externalFunctions.resize(extChunkSize);
    }

    void resize(const size_t& dataSize, const size_t& stmtSize) {
      data.resize(dataSize);
      statements.resize(stmtSize);
      adjoints.resize(stmtSize + 1);
    }

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
          assert(statements.getUsedSize() < statement.size);
          statements.setDataAndMove(std::make_tuple((StatementInt)activeVariables));
          lhsIndex = statements.getUsedSize();
        }
      }

      /* now set the value of the lhs */
      lhsValue = rhs.getValue();
    }

    inline void store(Real& value, IndexType& lhsIndex, const ActiveReal<Real, SimpleTape<Real, IndexType> >& rhs) {
      ENABLE_CHECK(OptTapeActivity, active){
        lhsIndex = rhs.getGradientData();
      }
      value = rhs.getValue();
    }

    inline void store(Real& value, IndexType& lhsIndex, const typename TypeTraits<Real>::PassiveReal& rhs) {
      ENABLE_CHECK(OptTapeActivity, active) {
        lhsIndex = 0;
      }
      value = rhs;
    }

    inline void pushJacobi(Real& /*gradient*/, const Real& /*value*/, const IndexType& index) {
      if(0 != index) {
        assert(data.getUsedSize() < data.size);

        data.setDataAndMove(std::make_tuple(1.0, index));
      }
    }

    inline void pushJacobi(Real& /*gradient*/, const Real& jacobi, const Real& /*value*/, const IndexType& index) {
      if(0 != index) {
        ENABLE_CHECK(OptJacobiIsZero, 0.0 != jacobi) {
          assert(data.getUsedSize() < data.size);

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
      assert(index < statements.size);
      return adjoints.data[index];
    }

    inline Real& gradient(IndexType& index) {
      assert(index < statements.size);
      assert(0 != index);

      return adjoints.data[index];
    }

    inline Position getPosition() {
      return Position(statements.getUsedSize(), data.getUsedSize(), externalFunctions.getUsedSize());
    }

    inline void reset(const Position& pos) {
      assert(pos.stmt < statements.size);
      assert(pos.data < data.size);

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

    inline void reset() {
      reset(Position(0,0,0));
    }

    inline void clearAdjoints(){
      for(size_t i = 0; i <= statements.getUsedSize(); ++i) {
        adjoints.data[i] = 0.0;
      }
    }

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

    inline void evaluate() {
      evaluate(getPosition(), Position(0,0,0));
    }

    inline void registerInput(ActiveReal<Real, SimpleTape<Real, IndexType> >& value) {
      assert(statements.getUsedSize() < statements.size);

      statements.setDataAndMove(std::make_tuple((StatementInt) 0));
      value.getGradientData() = statements.getUsedSize();
    }

    inline void registerOutput(ActiveReal<Real, SimpleTape<Real, IndexType> >& /*value*/) {
      /* do nothing */
    }

    inline void setActive(){
      active = true;
    }

    inline void setPassive(){
      active = false;
    }

    inline bool isActive(){
      ENABLE_CHECK(OptTapeActivity, true) {
        // default branch will return the tape activity
        return active;
      } else {
        // if we do not check for the tape activity, the tape is always active
        return true;
      }
    }

    void pushExternalFunctionHandle(ExternalFunction::CallFunction extFunc, void* checkpoint, ExternalFunction::DeleteFunction delCheckpoint){
      assert(0 != externalFunctions.getUnusedSize());
      ExternalFunction function(extFunc, checkpoint, delCheckpoint);
      externalFunctions.setDataAndMove(std::make_tuple(function, getPosition()));
    }

    template<typename Data>
    void pushExternalFunction(typename ExternalFunctionDataHelper<Data>::CallFunction extFunc, Data* checkpoint, typename ExternalFunctionDataHelper<Data>::DeleteFunction delCheckpoint){
      ExternalFunctionDataHelper<Data>* functionHelper = new ExternalFunctionDataHelper<Data>(extFunc, checkpoint, delCheckpoint);
      pushExternalFunctionHandle( ExternalFunctionDataHelper<Data>::callFunction, functionHelper, ExternalFunctionDataHelper<Data>::deleteFunction);
    }
  };
}


