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

#include "reuseIndexManagerBase.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Reuse index manager with a many-to-one relation between tapes and index manager.
   *
   * This is the classical implementation of the reuse index management strategy. For the details of reuse index
   * management, see ReuseIndexManagerBase.
   *
   * This index manager is not thread-safe.
   *
   * @tparam T_Index   Type for the identifier, usually an integer type.
   */
  template<typename T_Index>
  struct ReuseIndexManager : public ReuseIndexManagerBase<T_Index, ReuseIndexManager<T_Index>> {
    public:

      using Index = CODI_DD(T_Index, int);                           ///< See ReuseIndexManager.
      using Base = ReuseIndexManagerBase<Index, ReuseIndexManager>;  ///< Base class abbreviation.

      friend Base;  ///< Allow the base class to call protected and private methods.

      /*******************************************************************************/
      /// @name IndexManagerInterface: Constants
      /// @{

      using Base::CopyNeedsStatement;  ///< See ReuseIndexManagerBase.
      using Base::IsLinear;            ///< See ReuseIndexManagerBase.
      using Base::NeedsStaticStorage;  ///< See ReuseIndexManagerBase.

      /// @}

    private:

      Index globalMaximumIndex;  ///< The largest created index.

    public:

      /// Constructor
      ReuseIndexManager(Index const& reservedIndices) : globalMaximumIndex(reservedIndices) {
        generateNewIndices();
      }

      /// Destructor
      ~ReuseIndexManager() {}

      /*******************************************************************************/
      /// @name IndexManagerInterface: Methods
      /// @{

      /// \copydoc IndexManagerInterface::addToTapeValues <br><br>
      /// Implementation: Adds max live indices, cur live indices.
      void addToTapeValues(TapeValues& values) const {
        unsigned long maximumGlobalIndex = globalMaximumIndex;
        unsigned long storedIndices = this->usedIndicesPos + this->unusedIndicesPos;
        long currentLiveIndices = maximumGlobalIndex - storedIndices;

        values.addUnsignedLongEntry("Max. live indices", maximumGlobalIndex);
        values.addLongEntry("Cur. live indices", currentLiveIndices);

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

        for (size_t pos = 0; pos < this->indexSizeIncrement; ++pos) {
          this->unusedIndices[this->unusedIndicesPos + pos] = globalMaximumIndex + Index(pos) + 1;
        }

        this->unusedIndicesPos = this->indexSizeIncrement;
        globalMaximumIndex += this->indexSizeIncrement;
      }
  };
}
