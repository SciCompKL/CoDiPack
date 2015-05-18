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
    size_t op;
    size_t data;

    SimpleTapePosition(const size_t& op, const size_t& data) :
      op(op),
      data(data) {}
  };

  template <typename Real, typename IndexType>
  class SimpleTape : public ReverseTapeInterface<Real, IndexType, SimpleTape<Real, IndexType>, SimpleTapePosition > {
  public:

    typedef SimpleTapePosition Position;
  private:
    Chunk2<Real, IndexType> data;
    Chunk1<OperationInt> operators;
    Chunk1<Real> adjoints;

    bool active;

  public:
    SimpleTape() :
      data(0),
      operators(0),
      adjoints(1),
      active(false){}

    void resize(const size_t& dataSize, const size_t& opSize) {
      data.resize(dataSize);
      operators.resize(opSize);
      adjoints.resize(opSize + 1);
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
          assert(operators.getUsedSize() < operators.size);
          operators.setDataAndMove(std::make_tuple((OperationInt)activeVariables));
          lhsIndex = operators.getUsedSize();
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
      assert(index < operators.size);
      return adjoints.data[index];
    }

    inline Real& gradient(IndexType& index) {
      assert(index < operators.size);
      assert(0 != index);

      return adjoints.data[index];
    }

    inline Position getPosition() {
      return Position(operators.getUsedSize(), data.getUsedSize());
    }

    inline void reset(const Position& pos) {
      assert(pos.op < operators.size);
      assert(pos.data < data.size);

      for(size_t i = pos.op; i <= operators.getUsedSize(); ++i) {
        adjoints.data[i] = 0.0;
      }

      operators.setUsedSize(pos.op);
      data.setUsedSize(pos.data);
    }

    inline void reset() {
      reset(Position(0,0));
    }

    inline void clearAdjoints(){
      for(size_t i = 0; i <= operators.getUsedSize(); ++i) {
        adjoints.data[i] = 0.0;
      }
    }

    inline void evaluate(const Position& start, const Position& end) {
      assert(start.data >= end.data);
      assert(start.op >= end.op);

      Position curPos = start;

      while(curPos.op > end.op) {
        const Real& adj = adjoints.data[curPos.op];
        --curPos.op;
        const OperationInt& activeVariables = operators.data[curPos.op];
        ENABLE_CHECK(OptZeroAdjoint, adj != 0.0){
          for(OperationInt curVar = 0; curVar < activeVariables; ++curVar) {
            --curPos.data;

            adjoints.data[data.data2[curPos.data]] += adj * data.data1[curPos.data];
          }
        } else {
          curPos.data -= activeVariables;
        }
      }
    }

    inline void evaluate() {
      evaluate(getPosition(), Position(0,0));
    }

    inline void registerInput(ActiveReal<Real, SimpleTape<Real, IndexType> >& value) {
      assert(operators.getUsedSize() < operators.size);

      operators.setDataAndMove(std::make_tuple((OperationInt) 0));
      value.getGradientData() = operators.getUsedSize();
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

    /**
     * @brief Add a external function to the tape.
     *
     * The external function is called during the reverse evaluation of the tape. It can be used to
     * give special treatment to code sections which have simpler reverse implementation than the
     * AD tool.
     *
     * @param       extFunc The function which is called during the reverse evluation of the tape.
     * @param    checkpoint The data argument for the function. The tape takes procession of the data and will delete it.
     * @param delCheckpoint The delete function for the data.
     */
    void pushExternalFunctionHandle(ExternalFunction::CallFunction extFunc, void* checkpoint, ExternalFunction::DeleteFunction delCheckpoint) {}

    /**
     * @brief Add a external function to the tape.
     *
     * The external function is called during the reverse evaluation of the tape. It can be used to
     * give special treatment to code sections which have simpler reverse implementation than the
     * AD tool.
     *
     * @param       extFunc The function which is called during the reverse evluation of the tape.
     * @param    checkpoint The data argument for the function. The tape takes procession of the data and will delete it.
     * @param delCheckpoint The delete function for the data.
     *
     * @tparam Data The data type for the data.
     */
    template<typename Data>
    void pushExternalFunction(
        typename ExternalFunctionDataHelper<Data>::CallFunction extFunc,
        Data* checkpoint,
        typename ExternalFunctionDataHelper<Data>::DeleteFunction delCheckpoint) {}

  };
}


