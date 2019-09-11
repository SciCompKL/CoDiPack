/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2019 Chair for Scientific Computing (SciComp), TU Kaiserslautern
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
 * Authors: Max Sagebaum, Tim Albring, (SciComp, TU Kaiserslautern)
 */

#pragma once

#include "../../configure.h"
#include "misc.hpp"
#include <mutex>
#include <atomic>

/**
 * @brief Global namespace for CoDiPack - Code Differentiation Package
 */
namespace codi {
  template<typename IndexType>
  class ParallelGlobalIndexHandler {
    private:
      typedef IndexType Index;

      std::atomic<Index> nextIndex;

    public:
      ParallelGlobalIndexHandler() : nextIndex(1) {
      }

      CODI_INLINE Index getNextIndex() {
        return nextIndex.load();
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

      // creates an index chunk of default size
      CODI_INLINE void getRange(IndexRange<Index>& range) {
        getRange(getRangeSize(), range);
      }

      // creates an index chunk of custom size
      CODI_INLINE void getRange(const Index& size, IndexRange<Index>& range) {
        range.first = nextIndex.fetch_add(size);

        ENABLE_CHECK(IsOverflowCheck, range.first > range.first + size) {
          CODI_EXCEPTION("Overflow in global index handler. Use a larger index type or a reuse index manager.");
        }

        range.last = range.first + size - 1;
      }
  };
}


