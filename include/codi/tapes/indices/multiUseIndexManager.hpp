#pragma once

#include <vector>

#include "../../aux/macros.hpp"
#include "../../config.h"
#include "reuseIndexManager.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Extends the ReuseIndexManager with a copy optimization.
   *
   * Mathematical and implementational details are explained in \ref SBG2021Index.
   *
   * Performs reference counting for each index. If the reference count is zero, then the index is freed and given back
   * to the ReuseIndexManager.
   *
   * @tparam _Index   Type for the identifier, usually an integer type.
   */
  template<typename _Index>
  struct MultiUseIndexManager : public ReuseIndexManager<_Index> {
    public:

      using Index = CODI_DD(_Index, int);     ///< See MultiUseIndexManager.
      using Base = ReuseIndexManager<Index>;  ///< Base class abbreviation.

      /*******************************************************************************/
      /// @name IndexManagerInterface: Constants
      /// @{

      static bool constexpr CopyNeedsStatement =
          !Config::CopyOptimization;           ///< Copy optimization only active if configured.
      static bool constexpr IsLinear = false;  ///< Identifiers are not coupled to statements.

      /// @}

    private:

      std::vector<Index> indexUse;  ///< Reference counting for each index.

    public:

      /// Constructor
      MultiUseIndexManager(Index const& reservedIndices) : Base(reservedIndices), indexUse(Config::SmallChunkSize) {
        resizeUseVector();
      }

      /*******************************************************************************/
      /// @name ReuseIndexManager: Overwrites
      /// @{

      /// \copydoc ReuseIndexManager::addToTapeValues <br><br>
      /// Implementation: Additionally adds memory consumed by the index use vector.
      void addToTapeValues(TapeValues& values) const {
        Base::addToTapeValues(values);

        double memoryindexUseVector = (double)indexUse.size() * (double)(sizeof(Index));

        values.addDoubleEntry("Memory: index use vector", memoryindexUseVector, true, true);
      }

      /// \copydoc ReuseIndexManager::assignIndex
      CODI_INLINE bool assignIndex(Index& index) {
        bool generatedNewIndex = false;

        if (Base::InactiveIndex != index) {
          indexUse[index] -= 1;
        }

        if (Base::InactiveIndex != index && 0 == indexUse[index]) {
          indexUse[index] = 1;
          // index would be freed and used again so we keep it
        } else {
          index = Base::InactiveIndex;  // reset index here such that the base class will return a new one
          generatedNewIndex = Base::assignIndex(index);
          if (generatedNewIndex) {
            resizeUseVector();
          }

          indexUse[index] = 1;
        }

        return generatedNewIndex;
      }

      /// \copydoc ReuseIndexManager::assignUnusedIndex
      CODI_INLINE bool assignUnusedIndex(Index& index) {
        freeIndex(index);  // zero check is performed inside

        bool generatedNewIndex = Base::assignUnusedIndex(index);
        if (generatedNewIndex) {
          resizeUseVector();
        }

        indexUse[index] = 1;

        return generatedNewIndex;
      }

      /// \copydoc ReuseIndexManager::copyIndex
      CODI_INLINE void copyIndex(Index& lhs, Index const& rhs) {
        if (Config::CopyOptimization) {
          // Skip the logic if the indices are the same.
          // This also prevents the bug that if &lhs == &rhs, the left hand side will always be deactivated.
          if (lhs != rhs) {
            freeIndex(lhs);

            if (Base::InactiveIndex != rhs) {  // do not handle the zero index
              indexUse[rhs] += 1;

              lhs = rhs;
            }
          }
        } else {
          // path if copy optimizations are disabled
          assignIndex(lhs);
        }
      }

      /// \copydoc ReuseIndexManager::freeIndex
      CODI_INLINE void freeIndex(Index& index) {
        if (Base::valid && Base::InactiveIndex != index) {  // do not free the zero index
          indexUse[index] -= 1;

          if (indexUse[index] == 0) {  // only free the index if it is not used any longer
            Base::freeIndex(index);
          } else {
            index = Base::InactiveIndex;
          }
        }
      }

      /// @}

    private:

      CODI_NO_INLINE void resizeUseVector() {
        indexUse.resize(this->getLargestCreatedIndex() + 1);
      }
  };
}
