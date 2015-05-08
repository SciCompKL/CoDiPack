/**
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015 Chair for Scientific Computing, TU Kaiserslautern
 *
 * This file is part of CoDiPack.
 *
 * CoDiPack is free software: you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation, either version 2 of the
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
 * Authors: TODO
 */
#pragma once

#include <tuple>
#include <vector>

#include "chunk.hpp"

namespace codi {

  struct EmptyChunkVector {
    /**
     * @brief Position without any data.
     */
    struct Position {};

    inline Position getPosition() {
      return Position();
    }

    inline void reset(const Position& /*pos*/) {}
  };

  template<typename ChunkData, typename NestedVector = EmptyChunkVector>
  class ChunkVector {
  public:

    typedef typename NestedVector::Position NestedPosition;

    struct Position {
      size_t chunk;
      size_t data;

      NestedPosition inner;

      Position() :
        chunk(0),
        data(0),
        inner() {}

      Position(const size_t& chunk, const size_t& data, const NestedPosition& inner) :
        chunk(chunk),
        data(data),
        inner(inner) {}
    };

  private:
    std::vector<ChunkData* > chunks;
    std::vector<NestedPosition> positions;

    ChunkData* curChunk;
    size_t curChunkIndex;

    size_t chunkSize;

    NestedVector& nested;

  public:
    ChunkVector(const size_t& chunkSize, NestedVector& nested) :
      chunks(),
      positions(),
      curChunk(NULL),
      curChunkIndex(0),
      chunkSize(chunkSize),
      nested(nested)
    {
      curChunk = new ChunkData(chunkSize);
      chunks.push_back(curChunk);
      positions.push_back(nested.getPosition());
    }

    ~ChunkVector() {
      for(size_t i = 0; i < chunks.size(); ++i) {
        delete chunks[i];
      }
    }

    void setChunkSize(const size_t& chunkSize) {
      this->chunkSize = chunkSize;

      for(size_t i = 0; i < chunks.size(); ++i) {
        chunks[i]->resize(this->chunkSize);
      }
    }

    void resize(const size_t& totalSize) {
      size_t noOfChunks = totalSize / chunkSize;
      if(0 != totalSize % chunkSize) {
        noOfChunks += 1;
      }

      for(size_t i = chunks.size(); i < noOfChunks; ++i) {
        chunks.push_back(new ChunkData(chunkSize));
        positions.push_back(nested.getPosition());
      }
    }

    inline void nextChunk() {
      curChunk->store();

      curChunkIndex += 1;
      if(chunks.size() == curChunkIndex) {
        curChunk = new ChunkData(chunkSize);
        chunks.push_back(curChunk);
        positions.push_back(nested.getPosition());
      } else {
        curChunk = chunks[curChunkIndex];
        curChunk->reset();
        positions[curChunkIndex] = nested.getPosition();
      }
    }

    void reset(const Position& pos) {
      assert(pos.chunk < chunks.size());
      assert(pos.data < chunkSize);

      curChunk = chunks[pos.chunk];
      curChunk->load();
      curChunk->setUsedSize(pos.data);
      curChunkIndex = pos.chunk;

      nested.reset(pos.inner);
    }

    void reset() {
      reset(Position());
    }

    inline void reserveItems(const size_t items) {
      assert(items <= chunkSize);

      if(chunkSize < curChunk->getUsedSize() + items) {
        nextChunk();
      }
    }

    inline void setDataAndMove(const typename ChunkData::DataValues& data) {
      // this method should only be called if reserverItems has been called
      curChunk->setDataAndMove(data);
    }

    inline typename ChunkData::DataPointer getCurDataAtPos() {
      return curChunk->getDataPointer(curChunk->getUsedSize());
    }

    inline size_t getChunkPosition() {
      return curChunk->getUsedSize();
    }

    inline Position getPosition() {
      return Position(curChunkIndex, curChunk->getUsedSize(), nested.getPosition());
    }

    inline NestedPosition getInnerPosition(const size_t& chunkIndex) {
      return positions[chunkIndex];
    }

    inline typename ChunkData::DataPointer getDataAtPosition(const size_t& chunkIndex, const size_t& dataPos) {
      return chunks[chunkIndex]->dataPointer(dataPos);
    }

    inline size_t getChunkUsedData(const size_t& chunkIndex) {
      return chunks[chunkIndex]->getUsedSize();
    }
  };
}
