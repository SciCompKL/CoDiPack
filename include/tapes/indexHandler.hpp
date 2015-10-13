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

#include <vector>

namespace codi {

  template<typename IndexType>
  class IndexHandler {
    private:
      IndexType globalMaximumIndex;
      IndexType currentMaximumIndex;
      std::vector<IndexType> freeIndices;

    public:
      IndexHandler() :
        globalMaximumIndex(0),
        currentMaximumIndex(0),
        freeIndices() {}

      inline void freeIndex(IndexType& index) {
        if(0 != index) { // do not free the zero index
          if(currentMaximumIndex == index) {
            // freed index is the maximum one so we can decrease the count
            --currentMaximumIndex;
          } else {
            freeIndices.push_back(index);
          }

          index = 0;
        }
      }

      inline IndexType createIndex() {
        if(0 != freeIndices.size()) {
          IndexType index = freeIndices.back();
          freeIndices.pop_back();
          return index;
        } else {
          if(globalMaximumIndex == currentMaximumIndex) {
            ++globalMaximumIndex;
          }
          ++currentMaximumIndex;

          return currentMaximumIndex;
        }
      }

      inline void checkIndex(IndexType& index) {
        if(0 == index) {
          index = this->createIndex();
        }
      }

      inline void reset() {
        /* do nothing */
      }

      inline IndexType getMaximumGlobalIndex() {
        return globalMaximumIndex;
      }

      inline IndexType getCurrentIndex() {
        return globalMaximumIndex;
      }

      size_t getNumberStoredIndices() {
        return freeIndices.size();
      }

      size_t getNumberAllocatedIndices() {
        return freeIndices.capacity();
      }
  };
}
