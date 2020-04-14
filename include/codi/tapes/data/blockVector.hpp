#pragma once
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
  struct BlockVector : public DataInterface<_NestedVector> {
    public:

      using Chunk = DECLARE_DEFAULT(_Chunk, ChunkBase);
      using NestedVector = DECLARE_DEFAULT(_NestedVector, DataInterface);
      using InternalPosHandle = size_t;

      using NestedPosition = typename NestedVector::Position;

      using Position = ArrayPosition<NestedPosition>;

    private:
      Chunk chunk;

      NestedVector* nested;

    public:

      BlockVector(size_t const& chunkSize, NestedVector* nested) :
        chunk(chunkSize),
        nested(NULL)
      {
        setNested(nested);
      }

      BlockVector(size_t const& chunkSize) :
        chunk(chunkSize),
        nested(NULL)
      {}

      /*******************************************************************************
       * Section: Misc functions
       *
       * Description: TODO
       *
       */

      void addToTapeValues(TapeValues& values) const {
        size_t allocedSize    = chunk.getSize();
        size_t dataEntries    = getDataSize();
        size_t entrySize      = Chunk::EntrySize;

        double  memoryUsed  = (double)dataEntries*(double)entrySize* TapeValues::BYTE_TO_MB;
        double  memoryAlloc = (double)allocedSize*(double)entrySize* TapeValues::BYTE_TO_MB;

        values.addUnsignedLongEntry("Total number", dataEntries);
        values.addDoubleEntry("Memory used", memoryUsed, true, false);
        values.addDoubleEntry("Memory allocated", memoryAlloc, false, true);
      }

      template<typename TargetPosition>
      CODI_INLINE TargetPosition extractPosition(Position const& pos) const {
        return nested->template extractPosition<TargetPosition>(pos.inner);
      }

      CODI_INLINE Position extractPosition(Position const& pos) const {
        return pos;
      }

      CODI_INLINE size_t getDataSize() const {
        return chunk.getUsedSize();
      }

      CODI_INLINE Position getPosition() const {
        return Position(chunk.getUsedSize(), nested->getPosition());
      }

      CODI_INLINE size_t getPushedDataCount(InternalPosHandle const& startPos) {
        return chunk.getUsedSize() - startPos;
      }

      CODI_INLINE Position getZeroPosition() const {
        return Position(0, nested->getZeroPosition());
      }

      template<typename ... Data>
      CODI_INLINE void pushData(Data const& ... data) {
        // this method should only be called if reserveItems has been called
        chunk.pushData(data...);
      }

      CODI_INLINE InternalPosHandle reserveItems(size_t const& items) {
        codiAssert(chunk.getUsedSize() + items <= chunk.getSize());

        return chunk.getUsedSize();
      }

      void resize(size_t const& totalSize) {
        chunk.resize(totalSize);
      }

      void reset() {
        resetTo(getZeroPosition());
      }

      void resetHard() {
        chunk.resize(0);

        nested->resetHard();
      }

      void resetTo(Position const& pos) {
        codiAssert(pos.data <= chunk.getSize());

        chunk.setUsedSize(pos.data);

        nested->resetTo(pos.inner);
      }

      void setNested(NestedVector* v) {
        // Set nested is only called once during the initialization.
        codiAssert(NULL == this->nested);
        codiAssert(v->getZeroPosition() == v->getPosition());

        this->nested = v;
      }

      void swap(BlockVector<Chunk, NestedVector>& other) {
        chunk.swap(other.chunk);

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
        pHandle.setPointers(0, &chunk);

        size_t dataPos = start.data;
        pHandle.callNestedForward(
              /* arguments for callNestedForward */
              nested, dataPos, end.data,
              /* arguments for nested->evaluateForward */
              start.inner, end.inner, function, std::forward<Args>(args)...);

        codiAssert(dataPos == end.data);
      }

      template<typename Function, typename ... Args>
      CODI_INLINE void evaluateReverse(Position const& start, Position const& end,Function const& function,
                                       Args&&... args) {
        PointerStore<Chunk> pHandle;

        size_t dataPos = start.data;
        pHandle.setPointers(0, &chunk);

        pHandle.callNestedReverse(
              /* arguments for callNestedReverse */
              nested, dataPos, end.data,
              /* arguments for nested->evaluateReverse */
              start.inner, end.inner, function, std::forward<Args>(args)...);

        codiAssert(dataPos == end.data);
      }

      template<typename FunctionObject, typename ... Args>
      CODI_INLINE void forEachChunk(FunctionObject& function, bool recursive, Args&&... args) {

        function(&chunk, std::forward<Args>(args)...);

        if(recursive) {
          nested->forEachChunk(function, recursive, std::forward<Args>(args)...);
        }
      }

      template<typename FunctionObject, typename ... Args>
      CODI_INLINE void forEachForward(Position const& start, Position const& end, FunctionObject& function, Args&&... args) {
        codiAssert(start.data <= end.data);

        PointerStore<Chunk> pHandle;

        for(size_t dataPos = start.data; dataPos < end.data; dataPos += 1) {
          pHandle.setPointers(dataPos, &chunk);
          pHandle.call(function, std::forward<Args>(args)...);
        }
      }

      template<typename FunctionObject, typename ... Args>
      CODI_INLINE void forEachReverse(Position const& start, Position const& end, FunctionObject const& function, Args&&... args) {
        codiAssert(start.data >= end.data);

        PointerStore<Chunk> pHandle;

        // we do not initialize dataPos with start - 1 since the type can be unsigned
        for(size_t dataPos = start.data; dataPos > end.data; /* decrement is done inside the loop */) {
          dataPos -= 1; // decrement of loop variable

          pHandle.setPointers(dataPos, &chunk);
          pHandle.call(function, std::forward<Args>(args)...);
        }
      }
  };
}
