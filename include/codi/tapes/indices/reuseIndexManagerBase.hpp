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
#include "../data/emptyData.hpp"
#include "indexManagerInterface.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Identifiers are reused. Freed identifiers are assigned to new variables. Variables keep their indices as
   * long as they are active.
   *
   * This index manager does not implement a copy optimization. Therefore, every copy operation needs a statement,
   * but variables will keep their identifier as long as they are active.
   *
   * Mathematical and implementational details are explained in \ref SBG2021Index.
   *
   * For generalization reasons, it also extends from the EmptyData DataInterface.
   *
   * This class contains the basic logic for index reuse. The implementing class has to add a mechanism to generate new
   * indices.
   *
   * @tparam T_Index   Type for the identifier, usually an integer type.
   * @tparam T_Impl    Implementing class.
   */
  template<typename T_Index, typename T_Impl>
  struct ReuseIndexManagerBase : public IndexManagerInterface<T_Index>, public EmptyData {
    public:

      using Index = CODI_DD(T_Index, int);                ///< See ReuseIndexManagerBase.
      using Impl = CODI_DD(T_Impl, CODI_IMPLEMENTATION);  ///< See ReuseIndexManagerBase.
      using Base = IndexManagerInterface<Index>;          ///< Base class abbreviation.

      using Position = EmptyData::Position;  ///< See EmptyData.

      /*******************************************************************************/
      /// @name IndexManagerInterface: Constants
      /// @{

      static bool constexpr CopyNeedsStatement = true;  ///< No copy optimization is implemented.
      static bool constexpr IsLinear = false;           ///< Identifiers are not coupled to statements.
      static bool constexpr NeedsStaticStorage = true;  ///< Identifiers are managed globally.

      /// @}

    protected:

      std::vector<Index> usedIndices;  ///< Pool of indices that have already been used in this recording.
      size_t usedIndicesPos;           ///< Number of remaining used indices.

      std::vector<Index> unusedIndices;  ///< Pool of indices that have not been used in this recording yet.
      size_t unusedIndicesPos;           ///< Number of remaining unused indices.

      size_t indexSizeIncrement;  ///< Block size for index pool enlargement.

      bool valid;  ///< Prevent index free after destruction.

    private:
      /*******************************************************************************/
      /// @name General implementation
      /// @{

      /// Cast to the implementation.
      CODI_INLINE Impl& cast() {
        return static_cast<Impl&>(*this);
      }

      /// Method to generate new indices. Only called when unusedIndices is empty.
      CODI_NO_INLINE void generateNewIndices() {
        cast().generateNewIndices();
      }

      /// @}

    public:

      /// Constructor
      /// The constructor of the implementing class is expected to call generateNewIndices.
      ReuseIndexManagerBase()
          : usedIndices(),
            usedIndicesPos(0),
            unusedIndices(),
            unusedIndicesPos(0),
            indexSizeIncrement(Config::SmallChunkSize),
            valid(true) {
        increaseIndicesSize(unusedIndices);
      }

      /// Destructor
      ~ReuseIndexManagerBase() {
        valid = false;
      }

      /*******************************************************************************/
      /// @name IndexManagerInterface: Methods
      /// @{

      /// \copydoc codi::IndexManagerInterface::assignIndex
      template<typename Tape>
      CODI_INLINE bool assignIndex(Index& index) {
        bool generatedNewIndex = false;

        if (Base::InactiveIndex == index) {
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

        EventSystem<Tape>::notifyIndexAssignListeners(index);

        return generatedNewIndex;
      }

      /// \copydoc codi::IndexManagerInterface::assignUnusedIndex
      template<typename Tape>
      CODI_INLINE bool assignUnusedIndex(Index& index) {
        freeIndex<Tape>(index);  // Zero check is performed inside.

        bool generatedNewIndex = false;
        if (0 == unusedIndicesPos) {
          generateNewIndices();
          generatedNewIndex = true;
        }

        unusedIndicesPos -= 1;
        index = unusedIndices[unusedIndicesPos];

        EventSystem<Tape>::notifyIndexAssignListeners(index);

        return generatedNewIndex;
      }

      /// \copydoc codi::IndexManagerInterface::copyIndex
      template<typename Tape>
      CODI_INLINE void copyIndex(Index& lhs, Index const& rhs) {
        if (Base::InactiveIndex == rhs) {
          freeIndex<Tape>(lhs);
        } else {
          assignIndex<Tape>(lhs);
        }
      }

      /// \copydoc codi::IndexManagerInterface::freeIndex
      template<typename Tape>
      CODI_INLINE void freeIndex(Index& index) {
        if (valid && Base::InactiveIndex != index) {  // Do not free the zero index.

          EventSystem<Tape>::notifyIndexFreeListeners(index);

          if (usedIndicesPos == usedIndices.size()) {
            increaseIndicesSize(usedIndices);
          }

          usedIndices[usedIndicesPos] = index;
          usedIndicesPos += 1;

          index = Base::InactiveIndex;
        }
      }

      /// \copydoc codi::IndexManagerInterface::reset
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

      /// \copydoc codi::IndexManagerInterface::addToTapeValues <br><br>
      /// Implementation: Adds indices stored, memory used, memory allocated.
      void addToTapeValues(TapeValues& values) const {
        unsigned long storedIndices = this->usedIndicesPos + this->unusedIndicesPos;
        unsigned long allocatedIndices = this->usedIndices.size() + this->unusedIndices.size();

        double memoryStoredIndices = (double)storedIndices * (double)(sizeof(Index));
        double memoryAllocatedIndices = (double)allocatedIndices * (double)(sizeof(Index));

        values.addUnsignedLongEntry("Indices stored", storedIndices);
        values.addDoubleEntry("Memory used", memoryStoredIndices, true, false);
        values.addDoubleEntry("Memory allocated", memoryAllocatedIndices, false, true);
      }

      /// @}

    private:

      CODI_NO_INLINE void increaseIndicesSize(std::vector<Index>& v) {
        v.resize(v.size() + indexSizeIncrement);
      }

      CODI_NO_INLINE void increaseIndicesSizeTo(std::vector<Index>& v, size_t minimalSize) {
        codiAssert(v.size() < minimalSize);

        size_t increaseMul = (minimalSize - v.size()) / indexSizeIncrement + 1;  // +1 always rounds up.
        v.resize(v.size() + increaseMul * indexSizeIncrement);
      }
  };
}
