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

#include <algorithm>
#include <vector>

#include "../../config.h"
#include "../../misc/macros.hpp"
#include "../../tools/parallel/parallelToolbox.hpp"
#include "reuseIndexManagerBase.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Reuse index manager with a one-to-one relation between tapes and index manager.
   *
   * This is a thread-safe implementation of the reuse index management strategy. For the details of reuse index
   * management, see ReuseIndexManagerBase. The key difference is that multiple tape-local index managers can acquire
   * non-overlapping ranges of indices from the same global management.
   *
   * @tparam T_Index            Type for the identifier, usually an integer type.
   * @tparam T_ParallelToolbox  Tools used to make this index manager thread-safe.
   */
  template<typename T_Index, typename T_ParallelToolbox>
  struct ParallelReuseIndexManager
      : public ReuseIndexManagerBase<T_Index, ParallelReuseIndexManager<T_Index, T_ParallelToolbox>> {
    public:

      using Index = CODI_DD(T_Index, int);  ///< See ParallelReuseIndexManager.
      using ParallelToolbox = CODI_DD(T_ParallelToolbox,
                                      CODI_DEFAULT_PARALLEL_TOOLBOX);        ///< See ParallelReuseIndexManager.
      using Base = ReuseIndexManagerBase<Index, ParallelReuseIndexManager>;  ///< Base class abbreviation.
      friend Base;  ///< Allow the base class to access protected and private members.

    private:

      template<typename Type>
      using Atomic = typename ParallelToolbox::template Atomic<Type>;

      using ReadWriteMutex = typename ParallelToolbox::ReadWriteMutex;

    public:

      /*******************************************************************************/
      /// @name IndexManagerInterface: Constants
      /// @{

      using Base::CopyNeedsStatement;                    ///< See ReuseIndexManagerBase.
      using Base::IsLinear;                              ///< See ReuseIndexManagerBase.
      static bool constexpr NeedsStaticStorage = false;  ///< Indices are managed globally, but the global part is
                                                         ///< handled by static members of the manager itself.

      /// @}

    private:

      static Atomic<T_Index> globalMaximumIndex;      ///< The largest index created across all instances of this class.
      static bool globalMaximumIndexInitialized;      ///< Indicates whether globalMaximumIndex is initialized.
      static ReadWriteMutex globalMaximumIndexMutex;  ///< Safeguards globalMaximumIndex, globalMaximumIndexInitialized.

    public:

      /// Constructor
      /// For a tape class that uses this index manager, all tape instances are expected to pass the same number of
      /// reservedIndices to this constructor.
      ParallelReuseIndexManager(Index const& reservedIndices) {
        globalMaximumIndexMutex.lockWrite();
        if (!globalMaximumIndexInitialized) {
          globalMaximumIndex = reservedIndices;
          globalMaximumIndexInitialized = true;
        }
        globalMaximumIndexMutex.unlockWrite();
        generateNewIndices();
      }

      /// Destructor
      ~ParallelReuseIndexManager() {}

      /*******************************************************************************/
      /// @name IndexManagerInterface: Methods
      /// @{

      /// \copydoc IndexManagerInterface::addToTapeValues <br><br>
      /// Implementation: Adds max live indices, cur live indices, indices stored, memory used, memory allocated.
      void addToTapeValues(TapeValues& values) const {
        unsigned long maximumGlobalIndex = globalMaximumIndex;

        values.addUnsignedLongEntry("Max. live indices", maximumGlobalIndex);
        // The number of current live indices cannot be computed from one instance alone.
        // It equals the number of maximum live indices minus the number of indices stored across all instances.

        Base::addToTapeValues(values);
      }

      /// \copydoc IndexManagerInterface::getLargestCreatedIndex
      /// The following properties are specific to the ReuseIndexManager and inherited by the MultiUseIndexManager:
      /// 1. tape resets do not change the largest created index,
      /// 2. it is not guaranteed that the largest created index has been assigned to a variable already.
      CODI_INLINE Index getLargestCreatedIndex() const {
        return globalMaximumIndex;
      }

      /// @}

    private:

      CODI_NO_INLINE void generateNewIndices() {
        // This method is only called when unused indices are empty.
        // Initially, a number of unused indices is created which
        // equals the number of indices we generate now, therefore
        // we do not have to check for size.

        codiAssert(this->unusedIndices.size() >= this->indexSizeIncrement);

        Index upperIndexRangeBound = globalMaximumIndex += this->indexSizeIncrement;  // note: atomic operation
        Index lowerIndexRangeBound = upperIndexRangeBound - this->indexSizeIncrement;

        for (size_t pos = 0; pos < this->indexSizeIncrement; ++pos) {
          this->unusedIndices[this->unusedIndicesPos + pos] = lowerIndexRangeBound + Index(pos) + 1;
        }

        this->unusedIndicesPos = this->indexSizeIncrement;
      }
  };

  template<typename Index, typename ParallelToolbox>
  typename ParallelReuseIndexManager<Index, ParallelToolbox>::template Atomic<Index>
      ParallelReuseIndexManager<Index, ParallelToolbox>::globalMaximumIndex;

  template<typename Index, typename ParallelToolbox>
  bool ParallelReuseIndexManager<Index, ParallelToolbox>::globalMaximumIndexInitialized = false;

  template<typename Index, typename ParallelToolbox>
  typename ParallelReuseIndexManager<Index, ParallelToolbox>::ReadWriteMutex
      ParallelReuseIndexManager<Index, ParallelToolbox>::globalMaximumIndexMutex;
}
