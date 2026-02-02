/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2025 Chair for Scientific Computing (SciComp), University of Kaiserslautern-Landau
 * Homepage: http://scicomp.rptu.de
 * Contact:  Prof. Nicolas R. Gauger (codi@scicomp.uni-kl.de)
 *
 * Lead developers: Max Sagebaum, Johannes Blühdorn (SciComp, University of Kaiserslautern-Landau)
 *
 * This file is part of CoDiPack (http://scicomp.rptu.de/software/codi).
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

#include <cstdint>
#include <vector>

#include "../../config.h"
#include "../../misc/eventSystem.hpp"
#include "../../misc/macros.hpp"
#include "../../misc/mathUtility.hpp"
#include "../data/emptyData.hpp"
#include "indexManagerInterface.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /// Index tag pair struct stored in the active real.
  ///
  /// @tparam T_Index  The type of the identifier for the real value.
  /// @tparam T_Tag    The type of the tag value.
  template<typename T_Index, typename T_Tag>
  struct IndexTagPair {
      using Index = CODI_DD(T_Index, int);  ///< See IndexTagPair.
      using Tag = CODI_DD(T_Tag, int);      ///< See IndexTagPair.

      Index id;  ///< Identifier for the real value.
      Tag tag;   ///< Tag data for lifetime and life cycle detection.

      /// Equal comparison.
      CODI_INLINE bool operator==(IndexTagPair const& o) {
        return id == o.id && tag == o.tag;
      }

      /// Not equal comparison
      CODI_INLINE bool operator!=(IndexTagPair const& o) {
        return id != o.id || tag != o.tag;
      }
  };

  /**
   * @brief Implements the IndexManagerInterface and mimics a multi use index management.
   *
   * This index manager should not be used in production code.
   *
   * The index manager increases tag for each reset that is called. Therefore it can detect the usage
   * of old values.
   *
   * In addition, it remembers the use count of old values. It can therefore detect if an old value is released to
   * often.
   *
   * Otherwise it behaves like a linear index manager.
   *
   * @tparam T_Index   Type for the identifier, usually an integer type.
   */
  template<typename T_Index>
  struct DebugMultiUseIndexManager : public IndexManagerInterface<T_Index>, public EmptyData {
    public:

      using Index = CODI_DD(T_Index, int);  ///< See DebugMultiUseIndexManager.
      using Tag = std::uint8_t;             ///< Tag used for the recording live time management.

      using ActiveTypeIndexData = IndexTagPair<Index, Tag>;  ///< See IndexManagerInterface.
      using Base = IndexManagerInterface<T_Index>;           ///< Abbreviation for the base class.

    private:

      static size_t constexpr TAG_BITS = sizeof(std::uint8_t) * 8;
      static Index constexpr TAG_SIZE = 1 << TAG_BITS;
      static Index constexpr TAG_MASK = TAG_SIZE - 1;

    public:

      /*******************************************************************************/
      /// @name IndexManagerInterface: Constants
      /// @{

      static bool constexpr CopyNeedsStatement = !Config::CopyOptimization;  ///< See IndexManagerInterface.
      static bool constexpr IsLinear = false;                                ///< See IndexManagerInterface.
      static bool constexpr NeedsStaticStorage = true;                       ///< See IndexManagerInterface.

      /// @}

    private:

      bool valid;
      Index reservedIndices;
      Index nextNewIdentifier;
      std::uint8_t curTag = 1;                   // Do not use the zero tag.
      std::vector<std::vector<Index>> indexUse;  ///< Reference counting for each index.

    public:

      /// Constructor
      DebugMultiUseIndexManager(Index const& reservedIndices)
          : valid(true),
            reservedIndices(reservedIndices),
            nextNewIdentifier(reservedIndices + 1),
            indexUse(2, std::vector<Index>(Config::SmallChunkSize)) {
        resizeVectors();
      }

      /// Destructor
      ~DebugMultiUseIndexManager() {
        valid = false;
      }

      /*******************************************************************************/
      /// @name IndexManagerInterface: Interface implementation
      /// @{

      /// \copydoc IndexManagerInterface::addToTapeValues <br><br>
      /// Implementation: Adds size of index use vector.
      void addToTapeValues(TapeValues& values) const {
        double memoryindexUseVector = 0;
        for (std::vector<Index> const& cur : indexUse) {
          memoryindexUseVector += (double)indexUse.size() * (double)(sizeof(Index));
        }

        TapeValues::LocalReductionOperation constexpr operation =
            NeedsStaticStorage ? TapeValues::LocalReductionOperation::Max : TapeValues::LocalReductionOperation::Sum;

        values.addDoubleEntry("Memory: index use vector", memoryindexUseVector, operation, true, true);
      }

      /// \copydoc IndexManagerInterface::assignIndex
      template<typename Tape>
      CODI_INLINE bool assignIndex(ActiveTypeIndexData& data) {
        if (Base::InactiveIndex != data.id) {
          freeIndex<Tape>(data);
        }

        Index index = nextNewIdentifier;
        nextNewIdentifier += 1;

        if (getLargestCreatedIndex() < 0) {
          CODI_EXCEPTION("Overfow of identifiers.");
        }

        resizeVectors();

        data = {index, curTag};
        indexUse[curTag][index] = 1;

        return true;
      }

      /// \copydoc IndexManagerInterface::assignUnusedIndex
      template<typename Tape>
      CODI_INLINE bool assignUnusedIndex(ActiveTypeIndexData& index) {
        return assignIndex<Tape>(index);
      }

      /// \copydoc IndexManagerInterface::copyIndex
      template<typename Tape>
      CODI_INLINE void copyIndex(ActiveTypeIndexData& lhs, ActiveTypeIndexData const& rhs) {
        validateIdentifier(lhs);
        validateIdentifier(rhs, false);
        if (Config::CopyOptimization) {
          // Skip the logic if the indices are the same.
          // This also prevents the bug that if &lhs == &rhs, the left hand side will always be deactivated.
          if (lhs != rhs) {
            freeIndex<Tape>(lhs);

            if (Base::InactiveIndex != rhs.id) {  // Do not handle the zero index.
              EventSystem<Tape>::notifyIndexCopyListeners(rhs.id);

              indexUse[curTag][rhs.id] += 1;
              lhs = rhs;
            }
          }
        } else {
          // Path if copy optimizations are disabled.
          assignIndex<Tape>(lhs);
        }
      }

      /// \copydoc IndexManagerInterface::freeIndex
      template<typename Tape>
      CODI_INLINE void freeIndex(ActiveTypeIndexData& data) {
        validateIdentifier(data);

        if (valid && Base::InactiveIndex != data.id) {  // Do not free the zero index.
          indexUse[data.tag][data.id] -= 1;

          if (indexUse[data.tag][data.id] == 0) {  // Only free the index if it is not used any longer.
            EventSystem<Tape>::notifyIndexFreeListeners(data.id);
          }
        }

        data = {Base::InactiveIndex, 0};
      }

      /// \copydoc IndexManagerInterface::initIndex
      CODI_INLINE void initIndex(ActiveTypeIndexData& index) {
        index = ActiveTypeIndexData();
      }

      /// \copydoc codi::IndexManagerInterface::reset
      CODI_INLINE void reset() {
        nextNewIdentifier = 1 + reservedIndices;

        nextTag();
      }

      /// \copydoc codi::IndexManagerInterface::validateRhsIndex
      CODI_INLINE void validateRhsIndex(ActiveTypeIndexData const& data) const {
        validateIdentifier(data, false);
      }

      /// \copydoc IndexManagerInterface::getIndex
      CODI_INLINE Index const& getIndex(ActiveTypeIndexData const& data) {
        return data.id;
      }

      /// \copydoc IndexManagerInterface::getIndex
      CODI_INLINE Index& getIndex(ActiveTypeIndexData& data) {
        return data.id;
      }

      /// \copydoc codi::IndexManagerInterface::getLargestCreatedIndex
      CODI_INLINE Index getLargestCreatedIndex() const {
        return nextNewIdentifier - 1;
      }

      /// @}

    private:

      CODI_INLINE void nextTag() {
        // Switch to next tag
        curTag += 1;
        curTag = curTag % TAG_SIZE;

        if (curTag == 0) {
          curTag += 1;
        }

        // Clear append index use vector for the tag
        if (indexUse.size() <= curTag) {
          indexUse.push_back(std::vector<Index>(Config::SmallChunkSize));
        } else {
          indexUse[curTag].clear();
          resizeVectors();
        }
      }

      CODI_NO_INLINE void resizeVectors() {
        indexUse[curTag].resize(nextNewIdentifier);
      }

      CODI_INLINE void validateIdentifier(ActiveTypeIndexData const& data, bool isLhs = true) const {
        if (Base::InactiveIndex != data.id) {
          // Check validity of the tag.
          if (data.tag == 0 || data.tag >= indexUse.size()) {
            CODI_EXCEPTION("Invalid tag '%d' with index '%d'.", (int)data.tag, (int)data.id);
          }

          if (indexUse[data.tag][data.id] <= 0) {
            if (data.tag == curTag) {
              CODI_EXCEPTION("Index '%d(%d)' is used after it was finally deleted.", (int)data.id, (int)data.tag);
            } else {
              CODI_EXCEPTION("Deleted index '%d(%d)' from old iteration is used.", (int)data.id, (int)data.tag);
            }
          }

          if (data.tag != curTag) {
            if (!isLhs) {
              CODI_EXCEPTION("Index '%d' from an old iteration '%d' is used, current tag is %d.", (int)data.id,
                             (int)data.tag, (int)curTag);
            }
          }
        }
      }
  };
}
