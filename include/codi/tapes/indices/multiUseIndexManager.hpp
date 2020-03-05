#pragma once

#include <vector>

#include "../../aux/macros.h"
#include "../../config.h"
#include "reuseIndexManager.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  template<typename _Index>
  struct MultiUseIndexManager : public ReuseIndexManager<_Index> {
    public:

      using Index = DECLARE_DEFAULT(_Index, int);

      static bool const AssignNeedsStatement = !Config::AssignOptimization;
      static bool const IsLinear = false;

      using Base = ReuseIndexManager<Index>;

    private:

      std::vector<Index> indexUse;

    public:

      MultiUseIndexManager(Index const& reserveIndices) :
        Base(reserveIndices),
        indexUse(Config::SmallChunkSize)
      {
        resizeUseVector();
      }

      void addToTapeValues(TapeValues& values) const {

        Base::addToTapeValues(values);

        double memoryindexUseVector = (double)indexUse.size()*(double)(sizeof(Index)) * TapeValues::BYTE_TO_MB;

        values.addDoubleEntry("Memory index use vector", memoryindexUseVector, true, true);
      }

      CODI_INLINE bool assignIndex(Index& index) {
        bool generatedNewIndex = false;

        if(Base::UnusedIndex != index) {
          indexUse[index] -= 1;
        }

        if(Base::UnusedIndex != index && 0 == indexUse[index]) {
          indexUse[index] = 1;
          // index would be freed and used again so we keep it
        } else {

          index = Base::UnusedIndex; // Reset index here such that the base class will return a new one
          generatedNewIndex = Base::assignIndex(index);
          if(generatedNewIndex) {
            resizeUseVector();
          }

          indexUse[index] = 1;
        }

        return generatedNewIndex;
      }

      CODI_INLINE bool assignUnusedIndex(Index& index) {
        freeIndex(index); // zero check is performed inside

        bool generatedNewIndex = Base::assignUnusedIndex(index);
        if (generatedNewIndex) {
          resizeUseVector();
        }

        indexUse[index] = 1;

        return generatedNewIndex;
      }

      CODI_INLINE void copyIndex(Index& lhs, Index const& rhs) {
        if(Config::AssignOptimization) {
          // skip the logic if the indices are the same.
          // This also prevents the bug, that if &lhs == &rhs the left hand side will always be deactivated.
          if(lhs != rhs) {
            freeIndex(lhs);

            if(Base::UnusedIndex != rhs) { // do not handle the zero index
              indexUse[rhs] += 1;

              lhs = rhs;
            }
          }
        } else {
            // path if assign optimizations are disabled
            assignIndex(lhs);
        }
      }

      CODI_INLINE void freeIndex(Index& index) {
        if(Base::valid && Base::UnusedIndex != index) { // do not free the zero index
          indexUse[index] -= 1;

          if(indexUse[index] == 0) { // only free the index if it is not used any longer
            Base::freeIndex(index);
          } else {
            index = Base::UnusedIndex;
          }
        }
      }

    private:
      CODI_NO_INLINE void resizeUseVector() {
        indexUse.resize(this->getLargestAssignedIndex() + 1);
      }
  };
}
