#pragma once

#include <algorithm>
#include <vector>

#include "../../aux/macros.hpp"
#include "../../config.h"
#include "../data/emptyData.hpp"
#include "indexManagerInterface.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Identifiers are reused. Freed identifiers are given out to new variables. Variables with an identifier will
   * keep this one.
   *
   * This index manager does not implement a copy optimization. Therefore every copy operation needs a statement, but
   * variables will keep there identifier if they are always active.
   *
   * For generalization reasons it also extends from the EmptyData DataInterface.
   *
   * @tparam _Index   Type for the identifier usually an integer type.
   */
  template<typename _Index>
  struct ReuseIndexManager : public IndexManagerInterface<_Index>, public EmptyData {
    public:

      using Index = CODI_DD(_Index, int);  ///< See ReuseIndexManager
      using Base = IndexManagerInterface<Index>;  ///< Base class abbreviation

      /*******************************************************************************/
      /// @name IndexManagerInterface: Constants
      /// @{

      static bool constexpr CopyNeedsStatement = true;  ///< No copy optimization implemented.
      static bool constexpr IsLinear = false;  ///< Identifiers are not coupled to statements.

      /// @}

    private:

      Index globalMaximumIndex;

      std::vector<Index> usedIndices;
      size_t usedIndicesPos;

      std::vector<Index> unusedIndices;
      size_t unusedIndicesPos;

      size_t indexSizeIncrement;

    protected:

      bool valid;  ///< Prevent index free after destruction.

    public:

      /// Constructor
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

      /// Destructor
      ~ReuseIndexManager() {
        valid = false;
      }

      /*******************************************************************************/
      /// @name IndexManagerInterface: Methods
      /// @{

      /// \copydoc IndexManagerInterface::addToTapeValues <br><br>
      /// Implementation: Adds max live indices, cur live indices, indices stored, memory used, memory allocated.
      void addToTapeValues(TapeValues& values) const {

        unsigned long maximumGlobalIndex = globalMaximumIndex;
        unsigned long storedIndices      = usedIndicesPos + unusedIndicesPos;
        unsigned long allocatedIndices   = usedIndices.size() + unusedIndices.size();
        long currentLiveIndices          = maximumGlobalIndex - storedIndices;

        double memoryStoredIndices    = (double)storedIndices*(double)(sizeof(Index));
        double memoryAllocatedIndices = (double)allocatedIndices*(double)(sizeof(Index));

        values.addUnsignedLongEntry("Max. live indices", maximumGlobalIndex);
        values.addLongEntry("Cur. live indices", currentLiveIndices);
        values.addUnsignedLongEntry("Indices stored", storedIndices);
        values.addDoubleEntry("Memory used", memoryStoredIndices, true, false);
        values.addDoubleEntry("Memory allocated", memoryAllocatedIndices, false, true);
      }


      /// \copydoc IndexManagerInterface::assignIndex
      CODI_INLINE bool assignIndex(Index& index) {
        bool generatedNewIndex = false;

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

        return generatedNewIndex;
      }

      /// \copydoc IndexManagerInterface::assignUnusedIndex
      CODI_INLINE bool assignUnusedIndex(Index& index) {
        freeIndex(index); // zero check is performed inside

        bool generatedNewIndex = false;
        if (0 == unusedIndicesPos) {
          generateNewIndices();
          generatedNewIndex = true;
        }

        unusedIndicesPos -= 1;
        index = unusedIndices[unusedIndicesPos];

        return generatedNewIndex;
      }

      /// \copydoc IndexManagerInterface::copyIndex
      CODI_INLINE void copyIndex(Index& lhs, Index const& rhs) {
        if (Base::UnusedIndex == rhs) {
          freeIndex(lhs);
        } else {
          assignIndex(lhs);
        }
      }

      /// \copydoc IndexManagerInterface::freeIndex
      CODI_INLINE void freeIndex(Index& index) {
        if (valid && Base::UnusedIndex != index) { // do not free the zero index

          if (usedIndicesPos == usedIndices.size()) {
            increaseIndicesSize(usedIndices);
          }

          usedIndices[usedIndicesPos] = index;
          usedIndicesPos += 1;

          index = Base::UnusedIndex;
        }
      }

      /// \copydoc IndexManagerInterface::getLargestAssignedIndex
      CODI_INLINE Index getLargestAssignedIndex() const {
        return globalMaximumIndex;
      }

      /// \copydoc IndexManagerInterface::reset
      CODI_INLINE void reset() {
        size_t totalSize = usedIndicesPos + unusedIndicesPos;
        if (totalSize > unusedIndices.size()) {
          increaseIndicesSizeTo(unusedIndices, totalSize);
        }

        for (size_t pos = 0; pos < usedIndicesPos; ++pos) {
          unusedIndices[unusedIndicesPos + pos] = usedIndices[pos];
        }
        unusedIndicesPos = totalSize;
        usedIndicesPos = 0;

        if (Config::SortIndicesOnReset) {
          if (totalSize == unusedIndices.size()) {
            std::sort(unusedIndices.begin(), unusedIndices.end());
          } else {
            std::sort(&unusedIndices[0], &unusedIndices[unusedIndicesPos]);
          }
        }
      }

      /// @}
      /*******************************************************************************/
      /// @name EmptyData: Method overwrite
      /// @{

      /// \copydoc EmptyData::resetTo
      void resetTo(Position const& pos) {
        if (pos == getZeroPosition()) {
          reset();
        }
      }

      /// @}

    private:

      CODI_NO_INLINE void generateNewIndices() {
        // method is only called when unused indices are empty
        // initially it holds a number of unused indices which is
        // the same amount as we now generate, therefore we do not
        // check for size

        codiAssert(unusedIndices.size() >= indexSizeIncrement);

        for (size_t pos = 0; pos < indexSizeIncrement; ++pos) {
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
