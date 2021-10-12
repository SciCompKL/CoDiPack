/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2021 Chair for Scientific Computing (SciComp), TU Kaiserslautern
 * Homepage: http://www.scicomp.uni-kl.de
 * Contact:  Prof. Nicolas R. Gauger (codi@scicomp.uni-kl.de)
 *
 * Lead developers: Max Sagebaum, Johannes Blühdorn (SciComp, TU Kaiserslautern)
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
 *  - SciComp, TU Kaiserslautern:
 *    - Max Sagebaum
 *    - Johannes Blühdorn
 *    - Former members:
 *      - Tim Albring
 */
#pragma once

#include <algorithm>
#include <vector>

#include "../../misc/macros.hpp"
#include "../../tools/parallel/parallelToolbox.hpp"
#include "../../config.h"
#include "reuseIndexManagerBase.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Reuse index manager with a one-to-one relation between tapes and index manager.
   *
   * This is a parallel implementation of the reuse index management strategy. For the details of reuse index
   * management, see ReuseIndexManagerBase.
   *
   * This index manager is thread-safe.
   *
   * @tparam T_Index            Type for the identifier, usually an integer type.
   * @tparam T_ParallelToolbox  Tools used to make this index manager thread-safe.
   */
  template<typename T_Index, typename T_ParallelToolbox>
  struct ParallelReuseIndexManager : public ReuseIndexManagerBase<T_Index> {
    public:

      using Index = CODI_DD(T_Index, int);        ///< See ParallelReuseIndexManager.
      using ParallelToolbox = CODI_DD(T_ParallelToolbox, CODI_T(codi::ParallelToolbox<CODI_ANY, CODI_ANY>));  ///< See ParallelReuseIndexManager.
      using Base = IndexManagerInterface<Index>;  ///< Base class abbreviation.

    private:

      template<typename Type>
      using ParallelToolbox::template Atomic<Type>;
      using Mutex = typename ParallelToolbox::Mutex;
      using typename ParallelToolbox::Mutex;

    public:

      /*******************************************************************************/
      /// @name IndexManagerInterface: Constants
      /// @{

      using Base::CopyNeedsStatement;  ///< See ReuseIndexManagerBase.
      using Base::IsLinear;            ///< See ReuseIndexManagerBase.

      static bool constexpr IsThreadSafe = true;       ///< Is thread-safe.

      /// @}

    private:

      static Atomic<Index> globalMaximumIndex;    ///< The largest index created across all instances of this class.
      static bool globalMaximumIndexInitialized;  ///< Indicates whether globalMaximumIndex is initialized.
      static Mutex globalMaximumIndexMutex;       ///< Safeguards globalMaximumIndex, globalMaximumIndexInitialized.

    public:

      /// Constructor
      /// For a tape class that uses this index manager, all tape instances are expected to pass the same number of
      /// reservedIndices to this constructor.
      ParallelReuseIndexManager(Index const& reservedIndices) {
        globalMaximumIndexMutex.lock();
        if (!globalMaximumIndexInitialized) {
          globalMaximumIndex = reservedIndices;
          globalMaximumIndexInitialized = true;
        }
        globalMaximumIndexMutex.unlock();
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
        // number of current live indices cannot be computed from one instance alone
        // it equals maximum live indices minus number of indices stored across all instances

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

        Index upperIndexRangeBound = globalMaximumIndex += this->indexSizeIncrement; // note: atomic operation
        Index lowerIndexRangeBound = upperIndexRangeBound - this->indexSizeIncrement;

        for (size_t pos = 0; pos < this->indexSizeIncrement; ++pos) {
          this->unusedIndices[this->unusedIndicesPos + pos] = lowerIndexRangeBound + Index(pos) + 1;
        }

        this->unusedIndicesPos = this->indexSizeIncrement;
      }
  };

  template<typename Index, typename ParallelToolbox, typename Tape>
  ParallelToolbox::template Atomic<Index> ParallelReuseIndexManager<Index, ParallelToolbox, Tape>::globalMaximumIndex;

  template<typename Index, typename ParallelToolbox, typename Tape>
  bool ParallelReuseIndexManager<Index, ParallelToolbox, Tape>::globalMaximumIndexInitialized = false;

  template<typename Index, typename ParallelToolbox, typename Tape>
  typename ParallelToolbox::Mutex ParallelReuseIndexManager<Index, ParallelToolbox, Tape>::globalMaximumIndexMutex;
}
