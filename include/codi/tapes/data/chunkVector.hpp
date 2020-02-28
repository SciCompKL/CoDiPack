#pragma once

#include <vector>

#include "../../aux/macros.h"
#include "../../config.h"
#include "chunk.hpp"
#include "emptyVector.hpp"
#include "dataInterface.hpp"
#include "pointerStore.hpp"
#include "position.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  template<typename _Chunk, typename _NestedVector = EmptyVector>
  struct ChunkVector : public DataInterface<_NestedVector> {
    public:

      using Chunk = DECLARE_DEFAULT(_Chunk, ChunkBase);
      using NestedVector = DECLARE_DEFAULT(_NestedVector, DataInterface);

      using NestedPosition = typename NestedVector::Position;

      using Position = ChunkPosition<NestedPosition>;

    private:
      std::vector<Chunk*> chunks;
      std::vector<NestedPosition> positions;

      Chunk* curChunk;
      size_t curChunkIndex;

      size_t chunkSize;

      NestedVector* nested;

    public:

      ChunkVector(size_t const& chunkSize, NestedVector* nested) :
        chunks(),
        positions(),
        curChunk(NULL),
        curChunkIndex(0),
        chunkSize(chunkSize),
        nested(NULL)
      {
        setNested(nested);
      }

      ChunkVector(size_t const& chunkSize) :
        chunks(),
        positions(),
        curChunk(NULL),
        curChunkIndex(0),
        chunkSize(chunkSize),
        nested(NULL)
      {}

      ~ChunkVector() {
        for(size_t i = 0; i < chunks.size(); ++i) {
          delete chunks[i];
        }
      }

      /*******************************************************************************
       * Section: Misc functions
       *
       * Description: TODO
       *
       */

      template<typename TargetPosition>
      CODI_INLINE TargetPosition extractPosition(Position const& pos) const {
        return nested->template extractPosition<TargetPosition>(pos.inner);
      }

      CODI_INLINE Position extractPosition(Position const& pos) const {
        return pos;
      }

      CODI_INLINE size_t getDataSize() const {
        size_t size = 0;
        for(size_t i = 0; i < chunks.size(); ++i) {
          size += chunks[i]->getUsedSize();
        }

        return size;
      }

      CODI_INLINE Position getPosition() const {
        return Position(curChunkIndex, curChunk->getUsedSize(), nested->getPosition());
      }

      CODI_INLINE Position getZeroPosition() const {
        return Position(0, 0, nested->getZeroPosition());
      }

      template<typename ... Data>
      CODI_INLINE void pushData(Data const& ... data) {
        // this method should only be called if reserveItems has been called
        curChunk->pushData(data...);
      }

      CODI_INLINE void reserveItems(size_t const& items) {
        codiAssert(items <= chunkSize);

        if(chunkSize < curChunk->getUsedSize() + items) {
          nextChunk();
        }
      }

      void resize(size_t const& totalSize) {
        size_t noOfChunks = totalSize / chunkSize;
        if(0 != totalSize % chunkSize) {
          noOfChunks += 1;
        }

        for(size_t i = chunks.size(); i < noOfChunks; ++i) {
          chunks.push_back(new Chunk(chunkSize));
          positions.push_back(nested->getPosition());
        }
      }

      void reset() {
        resetTo(getZeroPosition());
      }

      void resetHard() {
        for(size_t i = 1; i < chunks.size(); ++i) {
          delete chunks[i];
        }

        chunks.resize(1);
        positions.resize(1);
        curChunk = chunks[0];

        curChunk->setUsedSize(0);
        curChunkIndex = 0;

        nested->resetHard();
      }

      void resetTo(Position const& pos) {
        codiAssert(pos.chunk < chunks.size());
        codiAssert(pos.data <= chunkSize);

        for(size_t i = pos.chunk + 1; i <= curChunkIndex; i += 1) {
          chunks[i]->reset();
        }

        curChunk = chunks[pos.chunk];
        curChunk->setUsedSize(pos.data);
        curChunkIndex = pos.chunk;

        nested->resetTo(pos.inner);
      }

      void setNested(NestedVector* v) {
        // Set nested is only called once during the initialization.
        codiAssert(NULL == this->nested);
        codiAssert(v->getZeroPosition() == v->getPosition());

        this->nested = v;

        curChunk = new Chunk(chunkSize);
        chunks.push_back(curChunk);
        positions.push_back(nested->getZeroPosition());
      }

      void swap(ChunkVector<Chunk, NestedVector>& other) {
        std::swap(chunks, other.chunks);
        std::swap(positions, other.positions);
        std::swap(curChunkIndex, other.curChunkIndex);
        std::swap(chunkSize, other.chunkSize);

        curChunk = chunks[curChunkIndex];
        other.curChunk = other.chunks[other.curChunkIndex];

        nested->swap(*other.nested);
      }

      /*******************************************************************************
       * Section: Iterator functions
       *
       * Description: TODO
       *
       */

      template<typename Function, typename ... Args>
      CODI_INLINE void evaluateForward(Position const& start, Position const& end,Function const& function,
                                       Args&&... args) {
        PointerStore<Chunk> pHandle;

        size_t dataPos = start.data;
        NestedPosition curInnerPos = start.inner;
        for(size_t curChunk = start.chunk; curChunk < end.chunk; curChunk += 1) {

          pHandle.setPointers(0, chunks[curChunk]);

          NestedPosition endInnerPos = positions[curChunk + 1];
          pHandle.callNestedForward(
                /* arguments for callNestedForward */
                nested, dataPos, chunks[curChunk]->getUsedSize(),
                /* arguments for nested->evaluateForward */
                curInnerPos, endInnerPos, function, std::forward<Args>(args)...);

          // After a full chunk is evaluated, the data position needs to be at the end of the chunk
          codiAssert(dataPos == chunks[curChunk]->getUsedSize());

          curInnerPos = endInnerPos;

          dataPos = 0;
        }

        // Iterate over the remainder also covers the case if the start chunk and end chunk are the same
        pHandle.setPointers(0, chunks[end.chunk]);
        pHandle.callNestedForward(
              /* arguments for callNestedForward */
              nested, dataPos, end.data,
              /* arguments for nested->evaluateReverse */
              curInnerPos, end.inner, function, std::forward<Args>(args)...);

        codiAssert(dataPos == end.data); // after the last chunk is evaluated, the data position needs to be at the end position
      }

      template<typename Function, typename ... Args>
      CODI_INLINE void evaluateReverse(Position const& start, Position const& end,Function const& function,
                                       Args&&... args) {
        PointerStore<Chunk> pHandle;

        size_t dataPos = start.data;
        NestedPosition curInnerPos = start.inner;
        for(size_t curChunk = start.chunk; curChunk > end.chunk; curChunk -= 1) {

          pHandle.setPointers(0, chunks[curChunk]);

          NestedPosition endInnerPos = positions[curChunk];
          pHandle.callNestedReverse(
                /* arguments for callNestedReverse */
                nested, dataPos, 0,
                /* arguments for nested->evaluateReverse */
                curInnerPos, endInnerPos, function, std::forward<Args>(args)...);

          codiAssert(dataPos == 0); // after a full chunk is evaluated, the data position needs to be zero

          curInnerPos = endInnerPos;

          dataPos = chunks[curChunk - 1]->getUsedSize();
        }

        // Iterate over the remainder also covers the case if the start chunk and end chunk are the same
        pHandle.setPointers(0, chunks[end.chunk]);
        pHandle.callNestedReverse(
              /* arguments for callNestedReverse */
              nested, dataPos, end.data,
              /* arguments for nested->evaluateReverse */
              curInnerPos, end.inner, function, std::forward<Args>(args)...);

        codiAssert(dataPos == end.data); // after the last chunk is evaluated, the data position needs to be at the end position
      }

      template<typename FunctionObject, typename ... Args>
      CODI_INLINE void forEachChunk(FunctionObject& function, bool recursive, Args&&... args) {

        for(size_t chunkPos = 0; chunkPos < chunks.size(); chunkPos += 1) {

          function(chunks[chunkPos], std::forward<Args>(args)...);
        }

        if(recursive) {
          nested->forEachChunk(function, recursive, std::forward<Args>(args)...);
        }
      }

      template<typename FunctionObject, typename ... Args>
      CODI_INLINE void forEachForward(Position const& start, Position const& end, FunctionObject& function, Args&&... args) {
        codiAssert(start.chunk < end.chunk || (start.chunk == end.chunk && start.data <= end.data));
        codiAssert(end.chunk < chunks.size());

        size_t dataStart = start.data;
        for(size_t chunkPos = start.chunk; chunkPos < end.chunk; chunkPos += 1) {

          forEachChunkEntryForward(chunkPos, dataStart, chunks[chunkPos]->getUsedSize(), function, std::forward<Args>(args)...);

          dataStart = 0;

        }

        forEachChunkEntryForward(end.chunk, dataStart, end.data, function, std::forward<Args>(args)...);
      }

      template<typename FunctionObject, typename ... Args>
      CODI_INLINE void forEachReverse(Position const& start, Position const& end, FunctionObject& function, Args&&... args) {
        codiAssert(start.chunk > end.chunk || (start.chunk == end.chunk && start.data >= end.data));
        codiAssert(start.chunk < chunks.size());

        size_t dataStart = start.data;
        for(size_t chunkPos = start.chunk; chunkPos > end.chunk; /* decrement is done inside the loop */) {

          forEachChunkEntryReverse(chunkPos, dataStart, 0, function, std::forward<Args>(args)...);

          chunkPos -= 1;
          dataStart = chunks[chunkPos]->getUsedSize(); // decrement of loop variable

        }

        forEachChunkEntryReverse(end.chunk, dataStart, end.data, function, std::forward<Args>(args)...);
      }

    private:

      template<typename FunctionObject, typename ... Args>
      CODI_INLINE void forEachChunkEntryForward(size_t const& chunkPos, size_t const& start, size_t const& end, FunctionObject& function, Args&&... args) {
        codiAssert(start <= end);
        codiAssert(chunkPos < chunks.size());

        PointerStore<Chunk> pHandle;

        for(size_t dataPos = start; dataPos < end; dataPos += 1) {
          pHandle.setPointers(dataPos, chunks[chunkPos]);
          pHandle.call(function, std::forward<Args>(args)...);
        }
      }

      template<typename FunctionObject, typename ... Args>
      CODI_INLINE void forEachChunkEntryReverse(size_t const& chunkPos, size_t const& start, size_t const& end,
                                          FunctionObject& function, Args&&... args) {
        codiAssert(start >= end);
        codiAssert(chunkPos < chunks.size());

        PointerStore<Chunk> pHandle;

        // we do not initialize dataPos with start - 1 since the type can be unsigned
        for(size_t dataPos = start; dataPos > end; /* decrement is done inside the loop */) {
          dataPos -= 1; // decrement of loop variable

          pHandle.setPointers(dataPos, chunks[chunkPos]);
          pHandle.call(function, std::forward<Args>(args)...);
        }
      }

      CODI_NO_INLINE void nextChunk() {
        curChunkIndex += 1;
        if(chunks.size() == curChunkIndex) {
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
}
