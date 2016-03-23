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

#include "../../configure.h"

/**
 * @brief Global namespace for CoDiPack - Code Differentiation Package
 */
namespace codi {

  /**
   * @brief Handles the indices that can be used and reused and counts the use of indices.
   *
   * New indices are generated if required. All the freed indices are stored in a list
   * and are reused when indices are needed.
   *
   * There is also a count for the indices that stores how often an index is used by a variable.
   * Therefor a tape that uses this handler does not need to write a statement if an assign operation
   * is evaluated.
   *
   * @tparam Index  The type for the handled indices.
   */
  template<typename Index>
  class ReuseIndexHandlerUseCount {
    public:
      /**
       * @brief The type definition for other tapes who want to access the type.
       */
      typedef Index IndexType;


      /**
       * @brief If it is required to write an assign statement after the index is copied.
       */
      const static bool AssignNeedsStatement = false;

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

      /**
       * @brief The vector that count for the indices how often they are used.
       */
      std::vector<Index> indexUse;

      /**
       * @brief The size increment for the index use vector.
       *
       * The size of the vector is incremented every time when the maximum global index does no longer fit in the vector.
       */
      size_t indexUseSizeIncrement;

    public:

      /**
       * @brief Create a handler that has no indices in use.
       */
      ReuseIndexHandlerUseCount() :
        globalMaximumIndex(0),
        currentMaximumIndex(0),
        freeIndices(),
        indexUseSizeIncrement(DefaultSmallChunkSize),
        indexUse(DefaultSmallChunkSize){}

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
          indexUse[index] -= 1;

          if(indexUse[index] == 0) { // only free the index if it not used any longer
            if(currentMaximumIndex == index) {
              // freed index is the maximum one so we can decrease the count
              --currentMaximumIndex;
            } else {
              freeIndices.push_back(index);
            }
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
        Index index;
        if(0 != freeIndices.size()) {
          index = freeIndices.back();
          freeIndices.pop_back();
        } else {
          if(globalMaximumIndex == currentMaximumIndex) {
            ++globalMaximumIndex;

            checkIndexUseSize();
          }
          index = ++currentMaximumIndex;
        }

        indexUse[index] = 1;

        return index;
      }

      /**
       * @brief Check if the index is active and if it is only used by this instance if not a new index is generated.
       *
       * @param[inout] index The current value of the index. If 0 then a new index is generated.
       */
      inline void assignIndex(Index& index) {
        if(0 == index) {
          index = this->createIndex();
        } else if(indexUse[index] > 1) {
          indexUse[index] -= 1;

          index = this->createIndex();
        }
      }

      inline void copyIndex(Index& lhs, const Index& rhs) {
        freeIndex(lhs);

        if(0 != rhs) { // do not handle the zero index
          indexUse[rhs] += 1;

          lhs = rhs;
        }
      }

      /**
       * @brief Not needed by this manager.
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
       *   the indices that are stored,
       *   the memory for the allocated indices and
       *   the memory for the index use vector.
       *
       * @param[in,out] out  The information is written to the stream.
       * @param[in]     hLine  The horizontal line that separates the sections of the output.
       *
       * @tparam Stream The type of the stream.
       */
      template<typename Stream>
      void printStatistics(Stream& out, const std::string hLine) const {
        size_t maximumGlobalIndex     = (size_t)this->getMaximumGlobalIndex();
        size_t storedIndices          = (size_t)this->getNumberStoredIndices();
        size_t currentLiveIndices     = (size_t)this->getCurrentIndex() - this->getNumberStoredIndices();

        double memoryStoredIndices    = (double)storedIndices*(double)(sizeof(Index)) * BYTE_TO_MB;
        double memoryIndexUse         = (double)this->indexUse.size()*(double)(sizeof(Index)) * BYTE_TO_MB;
        double memoryAllocatedIndices = (double)this->getNumberAllocatedIndices()*(double)(sizeof(Index)) * BYTE_TO_MB;

        out << hLine
            << "Indices\n"
            << hLine
            << "  Max. live indices:    " << std::setw(10) << maximumGlobalIndex << "\n"
            << "  Cur. live indices:    " << std::setw(10) << currentLiveIndices << "\n"
            << "  Indices stored:       " << std::setw(10) << storedIndices << "\n"
            << "  Memory used:          " << std::setiosflags(std::ios::fixed)
                                          << std::setprecision(2)
                                          << std::setw(10)
                                          << memoryStoredIndices << " MB" << "\n"
            << "  Memory allocated:     " << std::setiosflags(std::ios::fixed)
                                          << std::setprecision(2)
                                          << std::setw(10)
                                          << memoryAllocatedIndices << " MB" << "\n"
            << "  Memory index use vec: " << std::setiosflags(std::ios::fixed)
                                          << std::setprecision(2)
                                          << std::setw(10)
                                          << memoryIndexUse << " MB" << "\n";
      }

    private:

      /**
       * @brief Ensure that the index use vector is big enough for all indices.
       *
       * The method increases the size of the index use vector by the chunk
       * increment defined in the constructor.
       */
      inline void checkIndexUseSize() {
        if(indexUse.size() <= globalMaximumIndex) {
          this->indexUse.resize(indexUse.size() + indexUseSizeIncrement);
        }
      }
  };
}
