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

#include "activeReal.hpp"
#include "chunk.hpp"
#include "tapeInterface.hpp"

namespace codi {

  template <typename Real, typename IndexType>
  class SimpleTape : public TapeInterface<Real, IndexType> {
  public:

    struct Position {
      size_t op;
      size_t data;

      Position(const size_t& op, const size_t& data) :
        op(op),
        data(data) {}
    };

  private:
    Chunk2<Real, IndexType> data;
    Chunk1<IndexType> operators;
    Chunk1<Real> adjoints;

  public:
    SimpleTape() :
      data(0),
      operators(0),
      adjoints(1) {}

    void resize(const size_t& dataSize, const size_t& opSize) {
      data.resize(dataSize);
      operators.resize(opSize);
      adjoints.resize(opSize + 1);
    }

    template<typename Rhs>
    inline void store(Real& lhsValue, IndexType& lhsIndex, const Rhs& rhs) {
      Real gradient; /* This value will not be used */

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
        operators.data[operators.getUsedSize()] = activeVariables;
        lhsIndex = operators.increase();
      }

      /* now set the value of the lhs */
      lhsValue = rhs.getValue();
    }

    inline void store(Real& value, IndexType& lhsIndex, const ActiveReal<Real, SimpleTape<Real, IndexType> >& rhs) {
      lhsIndex = rhs.getGradientData();
      value = rhs.getValue();
    }

    inline void store(Real& value, IndexType& lhsIndex, const typename TypeTraits<Real>::PassiveReal& rhs) {
      lhsIndex = 0;
      value = rhs;
    }

    inline void pushJacobi(Real& /*gradient*/, const Real& /*value*/, const IndexType& index) {
      if(0 != index) {
        assert(data.getUsedSize() < data.size);

        data.data1[data.getUsedSize()] = 1.0;
        data.data2[data.getUsedSize()] = index;
        data.increase();
      }
    }

    inline void pushJacobi(Real& /*gradient*/, const Real& jacobi, const Real& /*value*/, const IndexType& index) {
      if(0 != index) {
        assert(data.getUsedSize() < data.size);

        data.data1[data.getUsedSize()] = jacobi;
        data.data2[data.getUsedSize()] = index;
        data.increase();
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
      assert(index <= operators.size);
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

    inline void evaluate(const Position& start, const Position& end) {
      assert(start.data >= end.data);
      assert(start.op >= end.op);

      Position curPos = start;

      while(curPos.op > end.op) {
        const Real& adj = adjoints.data[curPos.op];
        --curPos.op;
        const IndexType& activeVariables = operators.data[curPos.op];
        for(IndexType curVar = 0; curVar < activeVariables; ++curVar) {
          --curPos.data;

          adjoints.data[data.data2[curPos.data]] += adj * data.data1[curPos.data];
        }
      }
    }

    inline void evaluate() {
      evaluate(getPosition(), Position(0,0));
    }
    inline void registerInput(ActiveReal<Real, SimpleTape<Real, IndexType> >& value) {
      assert(operators.getUsedSize() < operators.size);

      operators.data[operators.getUsedSize()] = 0;
      value.getGradientData() = operators.increase();
    }

    inline void registerOutput(ActiveReal<Real, SimpleTape<Real, IndexType> >& /*value*/) {
      /* do nothing */
    }
  };
}


