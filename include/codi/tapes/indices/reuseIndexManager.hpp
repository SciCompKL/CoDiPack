#pragma once

#include <algorithm>
#include <vector>

#include "../../aux/macros.h"
#include "../../config.h"
#include "../data/emptyVector.hpp"
#include "indexManagerInterface.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  template<typename _Index>
  struct ReuseIndexManager : public IndexManagerInterface<_Index>, public EmptyVector {
    public:

      using Index = DECLARE_DEFAULT(_Index, int);
      using Base = IndexManagerInterface<Index>;

      static bool const AssignNeedsStatement = true;
      static bool const IsLinear = false;

    private:

      Index globalMaximumIndex;

      std::vector<Index> usedIndices;
      size_t usedIndicesPos;

      std::vector<Index> unusedIndices;
      size_t unusedIndicesPos;

      size_t indexSizeIncrement;

    protected:

      bool valid;

    public:

      ReuseIndexManager(Index const& reserveIndices) :
        globalMaximumIndex(reserveIndices + 1),
        usedIndices(),
        usedIndicesPos(0),
        unusedIndices(),
        unusedIndicesPos(0),
        indexSizeIncrement(Config::SmallChunkSize),
        valid(true)
      {
        increaseIndicesSize(unusedIndices);
        generateNewIndices();
      }

      ~ReuseIndexManager() {
        valid = false;
      }

      void addToTapeValues(TapeValues& values) const {

        unsigned long maximumGlobalIndex = globalMaximumIndex;
        unsigned long storedIndices      = usedIndicesPos + unusedIndicesPos;
        unsigned long allocatedIndices   = usedIndices.size() + unusedIndices.size();
        long currentLiveIndices          = maximumGlobalIndex - storedIndices;

        double memoryStoredIndices    = (double)storedIndices*(double)(sizeof(Index)) * TapeValues::BYTE_TO_MB;
        double memoryAllocatedIndices = (double)allocatedIndices*(double)(sizeof(Index)) * TapeValues::BYTE_TO_MB;

        values.addUnsignedLongEntry("Max. live indices", maximumGlobalIndex);
        values.addLongEntry("Cur. live indices", currentLiveIndices);
        values.addUnsignedLongEntry("Indices stored", storedIndices);
        values.addDoubleEntry("Memory used", memoryStoredIndices, true, false);
        values.addDoubleEntry("Memory allocated", memoryAllocatedIndices, false, true);
      }


      CODI_INLINE void assignIndex(Index& index, bool& generatedNewIndex = OptionalArg<bool>::value) {
        generatedNewIndex = false;

        if (Base::UnusedIndex == index) {
          if (0 == usedIndicesPos) {
            if (0 == unusedIndicesPos) {
              generateNewIndices();
              generatedNewIndex = true;
            }

            unusedIndicesPos -= 1;
            index = unusedIndices[unusedIndicesPos];
          } else {
            usedIndicesPos -= 1;
            index = usedIndices[usedIndicesPos];
          }
        }
      }

      CODI_INLINE void assignUnusedIndex(Index& index, bool& generatedNewIndex = OptionalArg<bool>::value) {
        freeIndex(index); // zero check is performed inside

        if(0 == unusedIndicesPos) {
          generateNewIndices();
          generatedNewIndex = true;
        } else {
          generatedNewIndex = false;
        }

        unusedIndicesPos -= 1;
        index = unusedIndices[unusedIndicesPos];
      }

      CODI_INLINE void copyIndex(Index& lhs, Index const& rhs) {
        if(Base::UnusedIndex == rhs) {
          freeIndex(lhs);
        } else {
          assignIndex(lhs);
        }
      }

      CODI_INLINE void freeIndex(Index& index) {
        if(valid && Base::UnusedIndex != index) { // do not free the zero index

          if(usedIndicesPos == usedIndices.size()) {
            increaseIndicesSize(usedIndices);
          }

          usedIndices[usedIndicesPos] = index;
          usedIndicesPos += 1;

          index = Base::UnusedIndex;
        }
      }

      CODI_INLINE Index getLargestAssignedIndex() const {
        return globalMaximumIndex;
      }

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

        if(Config::SortIndicesOnReset) {
          if(totalSize == unusedIndices.size()) {
            std::sort(unusedIndices.begin(), unusedIndices.end());
          } else {
            std::sort(&unusedIndices[0], &unusedIndices[unusedIndicesPos]);
          }
        }
      }

    private:

      CODI_NO_INLINE void generateNewIndices() {
        // method is only called when unused indices are empty
        // initially it holds a number of unused indices which is
        // the same amount as we now generate, therefore we do not
        // check for size

        codiAssert(unusedIndices.size() >= indexSizeIncrement);

        for(size_t pos = 0; pos < indexSizeIncrement; ++pos) {
          unusedIndices[unusedIndicesPos + pos] = globalMaximumIndex + Index(pos);
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
