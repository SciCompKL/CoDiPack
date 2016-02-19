/*
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

/**
 * @brief Global namespace for CoDiPack - Code Differentiation Package
 */
namespace codi {

  /**
   * @brief Handles the indices that can be used and reused.
   *
   * New indices are generated if required. All the freed indices are stored in a list
   * and are reused when indices are needed.
   *
   * @tparam Index  The type for the handled indices.
   */
  template<typename Index>
  class ReuseIndexHandler {
    public:
      /**
       * @brief The type definition for other tapes who want to access the type.
       */
      typedef Index IndexType;

    private:

      /** @brief The maximum index that was used over the whole process */
      Index globalMaximumIndex;
      /**
       * @brief The maximum index that is currently used.
       *
       * The maximum is decremented every time an index is released that corresponds to the
       * current maximum.
       */
      Index currentMaximumIndex;

      /**
       * @brief The list with the indices that are available for reuse.
       */
      std::vector<Index> freeIndices;

    public:

      /**
       * @brief Create a handler that has no indices in use.
       */
      ReuseIndexHandler() :
        globalMaximumIndex(0),
        currentMaximumIndex(0),
        freeIndices() {}

      /**
       * @brief Free the index that is given to the method.
       *
       * The method checks if the index is equal to the current maximum. If yes
       * then the current maximum is decremented otherwise the index is added
       * to the freed list.
       *
       * @param[inout] index  The index that is freed. It is set to zero in the method.
       */
      inline void freeIndex(Index& index) {
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

      /**
       * @brief Generate a new index.
       *
       * @return The new index that can be used.
       */
      inline Index createIndex() {
        if(0 != freeIndices.size()) {
          Index index = freeIndices.back();
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

      /**
       * @brief Check if the index is active if not a new index is generated.
       *
       * @param[inout] index The current value of the index. If 0 then a new index is generated.
       */
      inline void checkIndex(Index& index) {
        if(0 == index) {
          index = this->createIndex();
        }
      }

      /**
       * @brief Placeholder for further developments.
       */
      inline void reset() const {
        /* do nothing */
      }

      /**
       * @brief Get the maximum global
       *
       * @return The maximum index that was used during the lifetime of this index handler.
       */
      inline Index getMaximumGlobalIndex() const {
        return globalMaximumIndex;
      }

      /**
       * @brief Get the current maximum index.
       *
       * @return The current maximum index that is in use.
       */
      inline Index getCurrentIndex() const {
        return currentMaximumIndex;
      }

      /**
       * @brief Get the number of the stored indices.
       *
       * @return The number of stored indices.
       */
      size_t getNumberStoredIndices() const {
        return freeIndices.size();
      }

      /**
       * @brief Get the number of the allocated indices.
       *
       * @return The number of the allocated indices.
       */
      size_t getNumberAllocatedIndices() const {
        return freeIndices.capacity();
      }

      /**
       * @brief Output statistics about the used indices.
       *
       * Writes the 
       *   maximum number of live indices,
       *   the current number of lives indices,
       *   the indices that are stored and
       *   the memory for the allocated indices.
       */
      void printStatistics() const {
        size_t maximumGlobalIndex     = (size_t)this->getMaximumGlobalIndex();
        size_t storedIndices          = (size_t)this->getNumberStoredIndices();
        size_t currentLiveIndices     = (size_t)this->getCurrentIndex() - this->getNumberStoredIndices();

        double memoryStoredIndices    = (double)storedIndices*(double)(sizeof(Index)) * BYTE_TO_MB;
        double memoryAllocatedIndices = (double)this->getNumberAllocatedIndices()*(double)(sizeof(Index)) * BYTE_TO_MB;

        std::cout << "---------------------------------------------" << std::endl
                  << "Indices"                                       << std::endl
                  << "---------------------------------------------" << std::endl
                  << "  Max. live indices: " << std::setw(10) << maximumGlobalIndex << std::endl
                  << "  Cur. live indices: " << std::setw(10) << currentLiveIndices << std::endl
                  << "  Indices stored:    " << std::setw(10) << storedIndices << std::endl
                  << "  Memory allocated:  " << std::setiosflags(std::ios::fixed)
                                             << std::setprecision(2)
                                             << std::setw(10)
                                             << memoryAllocatedIndices << " MB" << std::endl
                  << "  Memory used:       " << std::setiosflags(std::ios::fixed)
                                             << std::setprecision(2)
                                             << std::setw(10)
                                             << memoryStoredIndices << " MB" << std::endl;
      }
  };
}
