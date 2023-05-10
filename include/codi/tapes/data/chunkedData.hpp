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
   * @brief Data is stored chunk-wise in this DataInterface implementation. If a chunk runs out of space, a new chunk
   * is allocated.
   *
   * See DataInterface documentation for details.
   *
   * Each chunk has the size provided in the constructor.
   *
   * @tparam T_Chunk            Has to implement ChunkBase. The chunk defines the data stored in this implementation.
   * @tparam T_NestedData       Nested DataInterface.
   * @tparam T_PointerInserter  Defines how data is appended to evaluate* function calls.
   */
  template<typename T_Chunk, typename T_NestedData = EmptyData, typename T_PointerInserter = PointerStore<T_Chunk>>
  struct ChunkedData : public DataInterface<T_NestedData> {
    public:

      using Chunk = CODI_DD(T_Chunk, Chunk1<CODI_ANY>);                                 ///< See ChunkedData
      using NestedData = CODI_DD(T_NestedData, CODI_T(DataInterface<CODI_ANY>));        ///< See ChunkedData
      using PointerInserter = CODI_DD(T_PointerInserter, CODI_T(PointerStore<Chunk>));  ///< See ChunkedData

      using InternalPosHandle = size_t;                      ///< Position in the chunk
      using NestedPosition = typename NestedData::Position;  ///< Position of NestedData

      using Position = ChunkPosition<NestedPosition>;  ///< \copydoc DataInterface::Position

    private:
      std::vector<Chunk*> chunks;
      std::vector<NestedPosition> positions;

      Chunk* curChunk;
      size_t curChunkIndex;

      size_t chunkSize;

      NestedData* nested;

    public:

      /// Allocate chunkSize entries and set the nested DataInterface.
      ChunkedData(size_t const& chunkSize, NestedData* nested)
          : chunks(), positions(), curChunk(nullptr), curChunkIndex(0), chunkSize(chunkSize), nested(nullptr) {
        setNested(nested);
      }

      /// Allocate chunkSize entries. Requires a call to #setNested.
      ChunkedData(size_t const& chunkSize)
          : chunks(), positions(), curChunk(nullptr), curChunkIndex(0), chunkSize(chunkSize), nested(nullptr) {}

      /// Destructor
      ~ChunkedData() {
        for (size_t i = 0; i < chunks.size(); ++i) {
          delete chunks[i];
        }
      }

      /*******************************************************************************/
      /// @name Adding items

      /// \copydoc DataInterface::pushData
      template<typename... Data>
      CODI_INLINE void pushData(Data const&... data) {
        // This method should only be called if reserveItems has been called.
        curChunk->pushData(data...);
      }

      /// \copydoc DataInterface::reserveItems <br><br>
      /// Implementation: Creates a new chunk if not enough space is left.
      CODI_INLINE InternalPosHandle reserveItems(size_t const& items) {
        codiAssert(items <= chunkSize);

        if (chunkSize < curChunk->getUsedSize() + items) {
          nextChunk();
        }

        return curChunk->getUsedSize();
      }

      /*******************************************************************************/
      /// @name Size management

      /// \copydoc DataInterface::resize
      void resize(size_t const& totalSize) {
        size_t noOfChunks = totalSize / chunkSize;
        if (0 != totalSize % chunkSize) {
          noOfChunks += 1;
        }

        for (size_t i = chunks.size(); i < noOfChunks; ++i) {
          chunks.push_back(new Chunk(chunkSize));
          positions.push_back(nested->getPosition());
        }
      }

      /// \copydoc DataInterface::reset
      void reset() {
        resetTo(getZeroPosition());
      }

      /// \copydoc DataInterface::resetHard
      void resetHard() {
        for (size_t i = 1; i < chunks.size(); ++i) {
          delete chunks[i];
        }

        chunks.resize(1);
        positions.resize(1);
        curChunk = chunks[0];

        curChunk->setUsedSize(0);
        curChunkIndex = 0;

        nested->resetHard();
      }

      /// \copydoc DataInterface::resetTo
      void resetTo(Position const& pos) {
        codiAssert(pos.chunk < chunks.size());
        codiAssert(pos.data <= chunkSize);

        for (size_t i = pos.chunk + 1; i <= curChunkIndex; i += 1) {
          chunks[i]->reset();
        }

        curChunk = chunks[pos.chunk];
        curChunk->setUsedSize(pos.data);
        curChunkIndex = pos.chunk;

        nested->resetTo(pos.inner);
      }

      /// \copydoc DataInterface::erase
      /// Implementation: If the given range start..end does not only overlap with parts of chunks but contains complete
      /// chunks, those completely contained chunks are deleted in the course of the erase.
      void erase(Position const& start, Position const& end, bool recursive = true) {
        size_t chunkRange = end.chunk - start.chunk;

        if (chunkRange == 0) {
          chunks[start.chunk]->erase(start.data, end.data);
        } else {
          // Treat first chunk.
          chunks[start.chunk]->erase(start.data, chunks[start.chunk]->getUsedSize());

          // Treat last chunk.
          chunks[end.chunk]->erase(0, end.data);

          // Erase completely covered chunks and free their memory. Covers also the case that there is no such chunk.
          chunks.erase(chunks.begin() + start.chunk + 1, chunks.begin() + end.chunk);
        }

        if (recursive) {
          nested->erase(start.inner, end.inner, recursive);
        }
      }

      /*******************************************************************************/
      /// @name Position functions

      /// \copydoc DataInterface::getDataSize
      CODI_INLINE size_t getDataSize() const {
        size_t size = 0;
        for (size_t i = 0; i < chunks.size(); ++i) {
          size += chunks[i]->getUsedSize();
        }

        return size;
      }

      /// \copydoc DataInterface::getPosition
      CODI_INLINE Position getPosition() const {
        return Position(curChunkIndex, curChunk->getUsedSize(), nested->getPosition());
      }

      /// \copydoc DataInterface::getPushedDataCount
      CODI_INLINE size_t getPushedDataCount(InternalPosHandle const& startPos) {
        return curChunk->getUsedSize() - startPos;
      }

      /// \copydoc DataInterface::getZeroPosition
      CODI_INLINE Position getZeroPosition() const {
        return Position(0, 0, nested->getZeroPosition());
      }

      /// \copydoc DataInterface::getDataPointers
      template<typename... Data>
      CODI_INLINE void getDataPointers(InternalPosHandle const& startPos, Data*&... data) {
        curChunk->dataPointer(startPos, data...);
      }

      /*******************************************************************************/
      /// @name Misc functions
      /// @{

      /// \copydoc DataInterface::addToTapeValues <br><br>
      /// Implementation: Adds: Total number, Number of chunks, Memory used, Memory allocated
      void addToTapeValues(TapeValues& values) const {
        size_t numberOfChunks = chunks.size();
        size_t dataEntries = getDataSize();
        size_t entrySize = Chunk::EntrySize;

        double memoryUsed = (double)dataEntries * (double)entrySize;
        double memoryAlloc = (double)numberOfChunks * (double)chunkSize * (double)entrySize;

        values.addUnsignedLongEntry("Total number", dataEntries);
        values.addUnsignedLongEntry("Number of chunks", numberOfChunks);
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

        curChunk = new Chunk(chunkSize);
        chunks.push_back(curChunk);
        positions.push_back(nested->getZeroPosition());
      }

      /// \copydoc DataInterface::swap
      void swap(ChunkedData<Chunk, NestedData>& other) {
        std::swap(chunks, other.chunks);
        std::swap(positions, other.positions);
        std::swap(curChunkIndex, other.curChunkIndex);
        std::swap(chunkSize, other.chunkSize);

        curChunk = chunks[curChunkIndex];
        other.curChunk = other.chunks[other.curChunkIndex];

        nested->swap(*other.nested);
      }

      /*******************************************************************************/
      /// @name Iterator functions

      /// \copydoc DataInterface::evaluateForward
      template<typename FunctionObject, typename... Args>
      CODI_INLINE void evaluateForward(Position const& start, Position const& end, FunctionObject function,
                                       Args&&... args) {
        PointerInserter pHandle;

        size_t curDataPos = start.data;
        size_t endDataPos;
        NestedPosition curInnerPos = start.inner;
        NestedPosition endInnerPos;

        size_t curChunk = start.chunk;
        for (;;) {
          // Update of end conditions.
          if (curChunk != end.chunk) {
            endInnerPos = positions[curChunk + 1];
            endDataPos = chunks[curChunk]->getUsedSize();
          } else {
            endInnerPos = end.inner;
            endDataPos = end.data;
          }

          pHandle.setPointers(0, chunks[curChunk]);
          pHandle.callNestedForward(
              /* arguments for callNestedForward */
              nested, curDataPos, endDataPos,
              /* arguments for nested->evaluateForward */
              curInnerPos, endInnerPos, function, std::forward<Args>(args)...);

          // After a full chunk is evaluated, the data position needs to be at the end data position.
          codiAssert(curDataPos == endDataPos);

          if (curChunk != end.chunk) {
            curChunk += 1;
            curInnerPos = endInnerPos;
            curDataPos = 0;
          } else {
            break;
          }
        }
      }

      /// \copydoc DataInterface::evaluateReverse
      template<typename FunctionObject, typename... Args>
      CODI_INLINE void evaluateReverse(Position const& start, Position const& end, FunctionObject function,
                                       Args&&... args) {
        PointerInserter pHandle;

        size_t curDataPos = start.data;
        size_t endDataPos;
        NestedPosition curInnerPos = start.inner;
        NestedPosition endInnerPos;

        size_t curChunk = start.chunk;
        for (;;) {
          // Update of end conditions.
          if (curChunk != end.chunk) {
            endInnerPos = positions[curChunk];
            endDataPos = 0;
          } else {
            endInnerPos = end.inner;
            endDataPos = end.data;
          }

          pHandle.setPointers(0, chunks[curChunk]);
          pHandle.callNestedReverse(
              /* arguments for callNestedReverse */
              nested, curDataPos, endDataPos,
              /* arguments for nested->evaluateReverse */
              curInnerPos, endInnerPos, function, std::forward<Args>(args)...);

          // After a full chunk is evaluated, the data position needs to be at the end data position.
          codiAssert(curDataPos == endDataPos);

          if (curChunk != end.chunk) {
            // Update of loop variables.
            curChunk -= 1;
            curInnerPos = endInnerPos;
            curDataPos = chunks[curChunk]->getUsedSize();
          } else {
            break;
          }
        }
      }

      /// \copydoc DataInterface::forEachChunk
      template<typename FunctionObject, typename... Args>
      CODI_INLINE void forEachChunk(FunctionObject& function, bool recursive, Args&&... args) {
        for (size_t chunkPos = 0; chunkPos < chunks.size(); chunkPos += 1) {
          function(chunks[chunkPos], std::forward<Args>(args)...);
        }

        if (recursive) {
          nested->forEachChunk(function, recursive, std::forward<Args>(args)...);
        }
      }

      /// \copydoc DataInterface::forEachForward
      template<typename FunctionObject, typename... Args>
      CODI_INLINE void forEachForward(Position const& start, Position const& end, FunctionObject function,
                                      Args&&... args) {
        codiAssert(start.chunk < end.chunk || (start.chunk == end.chunk && start.data <= end.data));
        codiAssert(end.chunk < chunks.size());

        size_t dataStart = start.data;
        for (size_t chunkPos = start.chunk; chunkPos <= end.chunk; chunkPos += 1) {
          size_t dataEnd;
          if (chunkPos != end.chunk) {
            dataEnd = chunks[chunkPos]->getUsedSize();
          } else {
            dataEnd = end.data;
          }

          forEachChunkEntryForward(chunkPos, dataStart, dataEnd, function, std::forward<Args>(args)...);

          dataStart = 0;
        }
      }

      /// \copydoc DataInterface::forEachReverse
      template<typename FunctionObject, typename... Args>
      CODI_INLINE void forEachReverse(Position const& start, Position const& end, FunctionObject function,
                                      Args&&... args) {
        codiAssert(start.chunk > end.chunk || (start.chunk == end.chunk && start.data >= end.data));
        codiAssert(start.chunk < chunks.size());

        size_t dataStart = start.data;
        size_t chunkPos = start.chunk;

        // For loop break condition is illformed due to unsigned underflow of chunkPos. The condition would be
        // chunkPos >= end.chunk which only breaks if chunkPos == -1 when end.chunk == 0. The minus one is not possible
        // for unsigned types.
        for (;;) {
          size_t dataEnd;
          if (chunkPos != end.chunk) {
            dataEnd = 0;
          } else {
            dataEnd = end.data;
          }

          forEachChunkEntryReverse(chunkPos, dataStart, dataEnd, function, std::forward<Args>(args)...);

          if (chunkPos == end.chunk) {
            break;
          } else {
            // Decrement of loop variable.
            chunkPos -= 1;
            dataStart = chunks[chunkPos]->getUsedSize();
          }
        }
      }

      /// @}

    private:

      template<typename FunctionObject, typename... Args>
      CODI_INLINE void forEachChunkEntryForward(size_t const& chunkPos, size_t const& start, size_t const& end,
                                                FunctionObject function, Args&&... args) {
        codiAssert(start <= end);
        codiAssert(chunkPos < chunks.size());

        PointerInserter pHandle;

        for (size_t dataPos = start; dataPos < end; dataPos += 1) {
          pHandle.setPointers(dataPos, chunks[chunkPos]);
          pHandle.call(function, std::forward<Args>(args)...);
        }
      }

      template<typename FunctionObject, typename... Args>
      CODI_INLINE void forEachChunkEntryReverse(size_t const& chunkPos, size_t const& start, size_t const& end,
                                                FunctionObject function, Args&&... args) {
        codiAssert(start >= end);
        codiAssert(chunkPos < chunks.size());

        PointerInserter pHandle;

        // For loop break condition is illformed due to unsigned underflow of dataPos. The condition would be
        // dataPos >= end which only breaks if dataPos == -1 when end == 0. The minus one is not possible
        // for unsigned types.
        for (size_t dataPos = start; dataPos > end; /* decrement is done inside the loop */) {
          dataPos -= 1;  // Decrement of loop variable.

          pHandle.setPointers(dataPos, chunks[chunkPos]);
          pHandle.call(function, std::forward<Args>(args)...);
        }
      }

      /// Loads next chunk or creates a new one if none is available.
      CODI_NO_INLINE void nextChunk() {
        curChunkIndex += 1;
        if (chunks.size() == curChunkIndex) {
          curChunk = new Chunk(chunkSize);
          chunks.push_back(curChunk);
          positions.push_back(nested->getPosition());
        } else {
          curChunk = chunks[curChunkIndex];
          curChunk->reset();
          positions[curChunkIndex] = nested->getPosition();
        }
      }
  };

  /// ChunkData DataInterface used in all regular tapes.
  template<typename Chunk, typename NestedData = EmptyData>
  using DefaultChunkedData = ChunkedData<Chunk, NestedData>;
}
