/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2023 Chair for Scientific Computing (SciComp), University of Kaiserslautern-Landau
 * Homepage: http://www.scicomp.uni-kl.de
 * Contact:  Prof. Nicolas R. Gauger (codi@scicomp.uni-kl.de)
 *
 * Lead developers: Max Sagebaum, Johannes Blühdorn (SciComp, University of Kaiserslautern-Landau)
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
 * For other licensing options please contact us.
 *
 * Authors:
 *  - SciComp, University of Kaiserslautern-Landau:
 *    - Max Sagebaum
 *    - Johannes Blühdorn
 *    - Former members:
 *      - Tim Albring
 */
#pragma once

#include <vector>

#include "../../config.h"
#include "../../misc/macros.hpp"
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
   * @tparam T_Index   Type for the identifier, usually an integer type.
   */
  template<typename T_Index>
  struct MultiUseIndexManager : public ReuseIndexManager<T_Index> {
    public:

      using Index = CODI_DD(T_Index, int);    ///< See MultiUseIndexManager.
      using Base = ReuseIndexManager<Index>;  ///< Base class abbreviation.

      /*******************************************************************************/
      /// @name IndexManagerInterface: Constants
      /// @{

      static bool constexpr CopyNeedsStatement =
          !Config::CopyOptimization;           ///< Copy optimization only active if configured.
      static bool constexpr IsLinear = false;  ///< See ReuseIndexManager.
      using Base::NeedsStaticStorage;          ///< See ReuseIndexManager.

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
      template<typename Tape>
      CODI_INLINE bool assignIndex(Index& index) {
        bool generatedNewIndex = false;

        if (Base::InactiveIndex != index) {
          indexUse[index] -= 1;
        }

        if (Base::InactiveIndex != index && 0 == indexUse[index]) {
          EventSystem<Tape>::notifyIndexFreeListeners(index);
          EventSystem<Tape>::notifyIndexAssignListeners(index);
          indexUse[index] = 1;
          // Index would be freed and used again so we keep it.
        } else {
          index = Base::InactiveIndex;  // Reset index here such that the base class will return a new one.
          generatedNewIndex = Base::template assignIndex<Tape>(index);
          if (generatedNewIndex) {
            resizeUseVector();
          }

          indexUse[index] = 1;
        }

        return generatedNewIndex;
      }

      /// \copydoc ReuseIndexManager::assignUnusedIndex
      template<typename Tape>
      CODI_INLINE bool assignUnusedIndex(Index& index) {
        freeIndex<Tape>(index);  // Zero check is performed inside.

        bool generatedNewIndex = Base::template assignUnusedIndex<Tape>(index);
        if (generatedNewIndex) {
          resizeUseVector();
        }

        indexUse[index] = 1;

        return generatedNewIndex;
      }

      /// \copydoc ReuseIndexManager::copyIndex
      template<typename Tape>
      CODI_INLINE void copyIndex(Index& lhs, Index const& rhs) {
        if (Config::CopyOptimization) {
          // Skip the logic if the indices are the same.
          // This also prevents the bug that if &lhs == &rhs, the left hand side will always be deactivated.
          if (lhs != rhs) {
            freeIndex<Tape>(lhs);

            if (Base::InactiveIndex != rhs) {  // Do not handle the zero index.
              EventSystem<Tape>::notifyIndexCopyListeners(rhs);

              indexUse[rhs] += 1;
              lhs = rhs;
            }
          }
        } else {
          // Path if copy optimizations are disabled.
          assignIndex<Tape>(lhs);
        }
      }

      /// \copydoc ReuseIndexManager::freeIndex
      template<typename Tape>
      CODI_INLINE void freeIndex(Index& index) {
        if (Base::valid && Base::InactiveIndex != index) {  // Do not free the zero index.
          indexUse[index] -= 1;

          if (indexUse[index] == 0) {  // Only free the index if it is not used any longer.
            Base::template freeIndex<Tape>(index);
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
