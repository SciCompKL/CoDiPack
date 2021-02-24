/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2021 Chair for Scientific Computing (SciComp), TU Kaiserslautern
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
 * Authors:
 *  - SciComp, TU Kaiserslautern:
 *     Max Sagebaum
 *     Tim Albring
 *     Johannes Blühdorn
 */

#pragma once

#include "../../configure.h"
#include "misc.hpp"

/**
 * @brief Global namespace for CoDiPack - Code Differentiation Package
 */
namespace codi {
  template<typename IndexType>
  class DefaultGlobalIndexHandler {
    private:
      typedef IndexType Index;

      Index nextIndex;

    public:
      DefaultGlobalIndexHandler() : nextIndex(1) {
      }

      CODI_INLINE Index getNextIndex() {
        return nextIndex;
      }

      CODI_INLINE Index getRangeSize() {
        return DefaultSmallChunkSize;
      }

      CODI_INLINE IndexRange<Index> getRange() {
        return getRange(getRangeSize());
      }

      CODI_INLINE IndexRange<Index> getRange(const Index& size) {
        IndexRange<Index> result;
        getRange(size, result);
        return result;
      }

      CODI_INLINE void getRange(IndexRange<Index>& range) {
        getRange(getRangeSize(), range);
      }

      CODI_INLINE void getRange(const Index& size, IndexRange<Index>& range) {
        ENABLE_CHECK(IsOverflowCheck, nextIndex > nextIndex + size) {
          CODI_EXCEPTION("Overflow in global index handler. Use a larger index type or a reuse index manager.");
        }

        range.first = nextIndex;
        range.last = nextIndex + size - 1;
        nextIndex += size;
      }
  };
}
