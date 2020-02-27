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

      CODI_INLINE void assignIndex(Index& index, bool& generatedNewIndex = OptionalArg<bool>::value) {
        generatedNewIndex = false;

        if(0 != index) {
          indexUse[index] -= 1;
        }

        if(0 != index && 0 == indexUse[index]) {
          indexUse[index] = 1;
          // index would be freed and used again so we keep it
        } else {

          index = 0; // Reset index here such that the base class will return a new one
          bool resized = false;
          Base::assignIndex(index, resized);
          if(resized) {
            resizeUseVector();
            generatedNewIndex = true;
          }

          indexUse[index] = 1;
        }
      }

      CODI_INLINE void assignUnusedIndex(Index& index, bool& generatedNewIndex = OptionalArg<bool>::value) {
        generatedNewIndex = false;

        freeIndex(index); // zero check is performed inside

        bool resized = false;
        Base::assignUnusedIndex(index, resized);

        if (resized) {
          resizeUseVector();
          generatedNewIndex = true;
        }

        indexUse[index] = 1;
      }

      CODI_INLINE void copyIndex(Index& lhs, Index const& rhs) {
        if(Config::AssignOptimization) {
          // skip the logic if the indices are the same.
          // This also prevents the bug, that if &lhs == &rhs the left hand side will always be deactivated.
          if(lhs != rhs) {
            freeIndex(lhs);

            if(0 != rhs) { // do not handle the zero index
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
        if(valid && 0 != index) { // do not free the zero index
          indexUse[index] -= 1;

          if(indexUse[index] == 0) { // only free the index if it is not used any longer
            Base::freeIndex(index);
          } else {
            index = 0;
          }
        }
      }

    private:
      CODI_NO_INLINE void resizeUseVector() {
        indexUse.resize(getLargestAssignedIndex() + 1);
      }
  };
}
