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
#include "../../traits/misc/enableIfHelpers.hpp"
#include "chunk.hpp"
#include "dataInterface.hpp"
#include "emptyData.hpp"
#include "pointerStore.hpp"
#include "position.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Data is stored in one contiguous block in this DataInterface implementation.
   *
   * See DataInterface documentation for details.
   *
   * This implementation does not check in #reserveItems if enough space is available. It needs to be preallocated with
   * #resize.
   *
   * @tparam T_Chunk            Has to implement ChunkBase. The chunk defines the data stored in this implementation.
   * @tparam T_NestedData       Nested DataInterface.
   * @tparam T_PointerInserter  Defines how data is appended to evaluate* function calls.
   */
  template<typename T_Chunk, typename T_NestedData = EmptyData, typename T_PointerInserter = PointerStore<T_Chunk>>
  struct BlockData : public DataInterface<T_NestedData> {
    public:

      using Chunk = CODI_DD(T_Chunk, CODI_T(Chunk1<CODI_ANY>));                         ///< See BlockData
      using NestedData = CODI_DD(T_NestedData, CODI_T(DataInterface<CODI_ANY>));        ///< See BlockData
      using PointerInserter = CODI_DD(T_PointerInserter, CODI_T(PointerStore<Chunk>));  ///< See BlockData

      using InternalPosHandle = size_t;                      ///< Position in the chunk
      using NestedPosition = typename NestedData::Position;  ///< Position of NestedData

      using Position = ArrayPosition<NestedPosition>;  ///< \copydoc DataInterface::Position

    private:
      Chunk chunk;

      NestedData* nested;

    public:

      /// Allocate chunkSize entries and set the nested DataInterface.
      BlockData(size_t const& chunkSize, NestedData* nested) : chunk(chunkSize), nested(nullptr) {
        setNested(nested);
      }

      /// Allocate chunkSize entries. Requires a call to #setNested.
      BlockData(size_t const& chunkSize) : chunk(chunkSize), nested(nullptr) {}

      /*******************************************************************************/
      /// @name Adding items

      /// \copydoc DataInterface::pushData
      template<typename... Data>
      CODI_INLINE void pushData(Data const&... data) {
        // This method should only be called if reserveItems has been called.
        chunk.pushData(data...);
      }

      /// \copydoc DataInterface::reserveItems <br><br>
      /// Implementation: Does not check if enough space is available.
      CODI_INLINE InternalPosHandle reserveItems(size_t const& items) {
        CODI_UNUSED(items);
        codiAssert(chunk.getUsedSize() + items <= chunk.getSize());

        return chunk.getUsedSize();
      }

      /*******************************************************************************/
      /// @name Size management

      /// \copydoc DataInterface::resize
      void resize(size_t const& totalSize) {
        chunk.resize(totalSize);
      }

      /// \copydoc DataInterface::reset
      void reset() {
        resetTo(getZeroPosition());
      }

      /// \copydoc DataInterface::resetHard
      void resetHard() {
        chunk.resize(0);

        nested->resetHard();
      }

      /// \copydoc DataInterface::resetTo
      void resetTo(Position const& pos) {
        codiAssert(pos.data <= chunk.getSize());

        chunk.setUsedSize(pos.data);

        nested->resetTo(pos.inner);
      }

      /// \copydoc DataInterface::erase
      /// Implementation: No memory is freed.
      void erase(Position const& start, Position const& end, bool recursive = true) {
        chunk.erase(start.data, end.data);
        if (recursive) {
          nested->erase(start.inner, end.inner);
        }
      }

      /*******************************************************************************/
      /// @name Position functions

      /// \copydoc DataInterface::getDataSize
      CODI_INLINE size_t getDataSize() const {
        return chunk.getUsedSize();
      }

      /// \copydoc DataInterface::getPosition
      CODI_INLINE Position getPosition() const {
        return Position(chunk.getUsedSize(), nested->getPosition());
      }

      /// \copydoc DataInterface::getPushedDataCount
      CODI_INLINE size_t getPushedDataCount(InternalPosHandle const& startPos) {
        return chunk.getUsedSize() - startPos;
      }

      /// \copydoc DataInterface::getZeroPosition
      CODI_INLINE Position getZeroPosition() const {
        return Position(0, nested->getZeroPosition());
      }

      /// \copydoc DataInterface::getDataPointers
      template<typename... Data>
      CODI_INLINE void getDataPointers(InternalPosHandle const& startPos, Data*&... data) {
        chunk.dataPointer(startPos, data...);
      }

      /*******************************************************************************/
      /// @name Misc functions

      /// \copydoc DataInterface::addToTapeValues <br><br>
      /// Implementation: Adds: Total number, Memory used, Memory allocated
      void addToTapeValues(TapeValues& values) const {
        size_t allocedSize = chunk.getSize();
        size_t dataEntries = getDataSize();
        size_t entrySize = Chunk::EntrySize;

        double memoryUsed = (double)dataEntries * (double)entrySize;
        double memoryAlloc = (double)allocedSize * (double)entrySize;

        values.addUnsignedLongEntry("Total number", dataEntries);
        values.addDoubleEntry("Memory used", memoryUsed, true, false);
        values.addDoubleEntry("Memory allocated", memoryAlloc, false, true);
      }

      /// \copydoc DataInterface::extractPosition
      template<typename TargetPosition, typename = typename enable_if_not_same<TargetPosition, Position>::type>
      CODI_INLINE TargetPosition extractPosition(Position const& pos) const {
        return nested->template extractPosition<TargetPosition>(pos.inner);
      }

      /// \copydoc DataInterface::extractPosition
      template<typename TargetPosition, typename = typename enable_if_same<TargetPosition, Position>::type>
      CODI_INLINE Position extractPosition(Position const& pos) const {
        return pos;
      }

      /// \copydoc DataInterface::setNested
      void setNested(NestedData* v) {
        // Set nested is only called once during the initialization.
        codiAssert(nullptr == this->nested);
        codiAssert(v->getZeroPosition() == v->getPosition());

        this->nested = v;
      }

      /// \copydoc DataInterface::swap
      void swap(BlockData<Chunk, NestedData>& other) {
        chunk.swap(other.chunk);

        nested->swap(*other.nested);
      }

      /*******************************************************************************/
      /// @name Iterator functions

      /// \copydoc DataInterface::evaluateForward
      template<typename FunctionObject, typename... Args>
      CODI_INLINE void evaluateForward(Position const& start, Position const& end, FunctionObject function,
                                       Args&&... args) {
        PointerInserter pHandle;
        pHandle.setPointers(0, &chunk);

        size_t dataPos = start.data;
        pHandle.callNestedForward(
            /* arguments for callNestedForward */
            nested, dataPos, end.data,
            /* arguments for nested->evaluateForward */
            start.inner, end.inner, function, std::forward<Args>(args)...);

        codiAssert(dataPos == end.data);
      }

      /// \copydoc DataInterface::evaluateReverse
      template<typename FunctionObject, typename... Args>
      CODI_INLINE void evaluateReverse(Position const& start, Position const& end, FunctionObject function,
                                       Args&&... args) {
        PointerInserter pHandle;

        size_t dataPos = start.data;
        pHandle.setPointers(0, &chunk);

        pHandle.callNestedReverse(
            /* arguments for callNestedReverse */
            nested, dataPos, end.data,
            /* arguments for nested->evaluateReverse */
            start.inner, end.inner, function, std::forward<Args>(args)...);

        codiAssert(dataPos == end.data);
      }

      /// \copydoc DataInterface::forEachChunk
      template<typename FunctionObject, typename... Args>
      CODI_INLINE void forEachChunk(FunctionObject& function, bool recursive, Args&&... args) {
        function(&chunk, std::forward<Args>(args)...);

        if (recursive) {
          nested->forEachChunk(function, recursive, std::forward<Args>(args)...);
        }
      }

      /// \copydoc DataInterface::forEachForward
      template<typename FunctionObject, typename... Args>
      CODI_INLINE void forEachForward(Position const& start, Position const& end, FunctionObject function,
                                      Args&&... args) {
        codiAssert(start.data <= end.data);

        PointerInserter pHandle;

        for (size_t dataPos = start.data; dataPos < end.data; dataPos += 1) {
          pHandle.setPointers(dataPos, &chunk);
          pHandle.call(function, std::forward<Args>(args)...);
        }
      }

      /// \copydoc DataInterface::forEachReverse
      template<typename FunctionObject, typename... Args>
      CODI_INLINE void forEachReverse(Position const& start, Position const& end, FunctionObject function,
                                      Args&&... args) {
        codiAssert(start.data >= end.data);

        PointerInserter pHandle;

        // We do not initialize dataPos with start - 1 since the type can be unsigned.
        for (size_t dataPos = start.data; dataPos > end.data; /* decrement is done inside the loop */) {
          dataPos -= 1;  // Decrement of loop variable.

          pHandle.setPointers(dataPos, &chunk);
          pHandle.call(function, std::forward<Args>(args)...);
        }
      }
  };

  /// BlockData DataInterface used in all unchecked tapes.
  template<typename Chunk, typename NestedData = EmptyData>
  using DefaultBlockData = BlockData<Chunk, NestedData>;
}
