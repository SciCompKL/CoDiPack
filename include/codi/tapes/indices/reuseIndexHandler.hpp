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

#include <vector>

#include "../../configure.h"
#include "../../tools/tapeValues.hpp"

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
   * @tparam IndexType  The type for the handled indices.
   */
  template<typename IndexType>
  class ReuseIndexHandler {
    public:
      /**
       * @brief The type definition for other tapes who want to access the type.
       */
      typedef IndexType Index;


      /**
       * @brief If it is required to write an assign statement after the index is copied.
       */
      const static bool AssignNeedsStatement = true;

      /**
       * @brief Indicates if the index handler privides linear increasing indices.
       *
       * false for this index manager.
       */
      static const bool IsLinear = false;

    private:

      /** @brief The maximum index that was used over the whole process */
      Index globalMaximumIndex;

      /**
       * @brief The list with the indices that are available for reuse.
       *
       * These indices have already been used.
       */
      std::vector<Index> usedIndices;

      /**
       * @brief The position for the used indices.
       */
      size_t usedIndicesPos;

      /**
       * @brief The list with the indices that are available for reuse.
       *
       * These indices have not been used in the current tape.
       */
      std::vector<Index> unusedIndices;

      /**
       * @brief The position for the unused indices.
       */
      size_t unusedIndicesPos;

      /**
       * @brief The number of indices generated if new indices are required.
       */
      size_t indexSizeIncrement;

      /**
       * @brief Indicates the destruction of the index handler.
       *
       * Required to prevent segmentation faults if varaibles are deleted after the index handler.
       */
      bool valid;

    public:

      /**
       * @brief Create a handler for index reuse.
       *
       * The argument reserveIndices will cause the index manager to reserve the first n indices, so that there are
       * not used by the index manager and are freely available to anybody.
       *
       * @param[in] reserveIndices  The number of indices that are reserved and not used by the manager.
       */
      ReuseIndexHandler(const Index reserveIndices) :
        globalMaximumIndex(reserveIndices + 1),
        usedIndices(),
        usedIndicesPos(0),
        unusedIndices(),
        unusedIndicesPos(0),
        indexSizeIncrement(DefaultSmallChunkSize),
        valid(true)
      {
        increaseIndicesSize(unusedIndices);
        generateNewIndices();
      }

      ~ReuseIndexHandler() {
        valid = false;
      }

      /**
       * @brief Free the index that is given to the method.
       *
       * The method checks if the index is equal to the current maximum. If yes
       * then the current maximum is decremented otherwise the index is added
       * to the freed list.
       *
       * @param[in,out] index  The index that is freed. It is set to zero in the method.
       */
      CODI_INLINE void freeIndex(Index& index) {
        if(valid && 0 != index) { // do not free the zero index

#if CODI_IndexHandle
          handleIndexFree(index);
#endif
          if(usedIndicesPos == usedIndices.size()) {
            increaseIndicesSize(usedIndices);
          }

          usedIndices[usedIndicesPos] = index;
          usedIndicesPos += 1;

          index = 0;
        }
      }

      /**
       * @brief Generate a new index.
       *
       * @return The new index that can be used.
       */
      CODI_INLINE Index createIndex() {
        Index index;

        if(0 == usedIndicesPos) {
          if(0 == unusedIndicesPos) {
            generateNewIndices();
          }

          unusedIndicesPos -= 1;
          index = unusedIndices[unusedIndicesPos];
        } else {
          usedIndicesPos -= 1;
          index = usedIndices[usedIndicesPos];
        }

#if CODI_IndexHandle
        handleIndexCreate(index);
#endif

        return index;
      }

      /**
       * @brief Generate a new index that has not been used in the tape.
       *
       * @return The new index that can be used.
       */
      CODI_INLINE Index createUnusedIndex() {
        Index index;

        if(0 == unusedIndicesPos) {
          generateNewIndices();
        }

        unusedIndicesPos -= 1;
        index = unusedIndices[unusedIndicesPos];

#if CODI_IndexHandle
        handleIndexCreate(index);
#endif

        return index;
      }

      /**
       * @brief Check if the index is active if not a new index is generated.
       *
       * @param[in,out] index The current value of the index. If 0 then a new index is generated.
       */
      CODI_INLINE void assignIndex(Index& index) {
        if(0 == index) {
          index = this->createIndex();
        }
      }

      /**
       * @brief Check if the index is active if yes it is deleted and an unused index is generated.
       *
       * @param[in,out] index The current value of the index.
       */
      CODI_INLINE void assignUnusedIndex(Index& index) {
        freeIndex(index); // zero check is performed inside

        index = this->createUnusedIndex();
      }

      /**
       * @brief The index on the rhs is ignored in this manager.
       *
       * The manager ensures only that the lhs index is valid.
       */
      CODI_INLINE void copyIndex(Index& lhs, const Index& rhs) {
        CODI_UNUSED(rhs);
        assignIndex(lhs);
      }

      /**
       * @brief Adds all used indices to the unused indices vector.
       */
      CODI_INLINE void reset() {
        size_t totalSize = usedIndicesPos + unusedIndicesPos;
        if(totalSize > unusedIndices.size()) {
          increaseIndicesSizeTo(unusedIndices, totalSize);
        }

        for(size_t pos = 0; pos < usedIndicesPos; ++pos) {
          unusedIndices[unusedIndicesPos + pos] = usedIndices[pos];
        }
        unusedIndicesPos = totalSize;
        usedIndicesPos = 0;

        if(OptSortIndicesOnReset) {
          if(totalSize == unusedIndices.size()) {
            std::sort(unusedIndices.begin(), unusedIndices.end());
          } else {
            std::sort(&unusedIndices[0], &unusedIndices[unusedIndicesPos]);
          }
        }
      }

      /**
       * @brief Get the maximum global
       *
       * @return The maximum index that was used during the lifetime of this index handler.
       */
      CODI_INLINE Index getMaximumGlobalIndex() const {
        return globalMaximumIndex;
      }

      /**
       * @brief Get the current maximum index.
       *
       * @return The current maximum index that is in use.
       */
      CODI_INLINE Index getCurrentIndex() const {
        return globalMaximumIndex;
      }

      /**
       * @brief Get the number of the stored indices.
       *
       * @return The number of stored indices.
       */
      CODI_INLINE size_t getNumberStoredIndices() const {
        return unusedIndicesPos + usedIndicesPos;
      }

      /**
       * @brief Get the number of the allocated indices.
       *
       * @return The number of the allocated indices.
       */
      CODI_INLINE size_t getNumberAllocatedIndices() const {
        return unusedIndices.capacity() + usedIndices.capacity();
      }

      /**
       * @brief Add statistics about the used indices.
       *
       * Adds the
       *   maximum number of live indices,
       *   the current number of lives indices,
       *   the indices that are stored and
       *   the memory for the allocated indices.
       *
       * @param[in,out] values  The values where the information is added to.
       */
      void addValues(TapeValues& values) const {
        size_t maximumGlobalIndex     = (size_t)this->getMaximumGlobalIndex();
        size_t storedIndices          = (size_t)this->getNumberStoredIndices();
        size_t currentLiveIndices     = (size_t)this->getCurrentIndex() - this->getNumberStoredIndices();

        double memoryStoredIndices    = (double)storedIndices*(double)(sizeof(Index)) * BYTE_TO_MB;
        double memoryAllocatedIndices = (double)this->getNumberAllocatedIndices()*(double)(sizeof(Index)) * BYTE_TO_MB;

        values.addSection("Indices");
        values.addData("Max. live indices", maximumGlobalIndex);
        values.addData("Cur. live indices", currentLiveIndices);
        values.addData("Indices stored", storedIndices);
        values.addData("Memory used", memoryStoredIndices, true, false);
        values.addData("Memory allocated", memoryAllocatedIndices, false, true);
      }

    private:
      CODI_NO_INLINE void generateNewIndices() {
        // method is only called when unused indices are empty
        // initialy it holds a number of unused indices which is
        // the same amount as we now generate, therefore we do not
        // check for size

        codiAssert(unusedIndices.size() >= indexSizeIncrement);

        for(size_t pos = 0; pos < indexSizeIncrement; ++pos) {
          unusedIndices[unusedIndicesPos + pos] = globalMaximumIndex + (Index)pos;
        }

        unusedIndicesPos = indexSizeIncrement;
        globalMaximumIndex += indexSizeIncrement;
      }

      CODI_NO_INLINE void increaseIndicesSize(std::vector<Index>& v) {
        v.resize(v.size() + indexSizeIncrement);
      }

      CODI_NO_INLINE void increaseIndicesSizeTo(std::vector<Index>& v, size_t minimalSize) {
        codiAssert(v.size() < minimalSize);

        size_t increaseMul = (minimalSize - v.size()) / indexSizeIncrement + 1; // +1 rounds always up
        v.resize(v.size() + increaseMul * indexSizeIncrement);
      }
  };
}
