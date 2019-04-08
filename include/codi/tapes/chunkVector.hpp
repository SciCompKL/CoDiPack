/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2019 Chair for Scientific Computing (SciComp), TU Kaiserslautern
 * Homepage: http://www.scicomp.uni-kl.de
 * Contact:  Prof. Nicolas R. Gauger (codi@scicomp.uni-kl.de)
 *
 * Lead developers: Max Sagebaum, Tim Albring (SciComp, TU Kaiserslautern)
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
 * Authors: Max Sagebaum, Tim Albring, (SciComp, TU Kaiserslautern)
 */

#pragma once

#include <tuple>
#include <vector>

#include "chunk.hpp"
#include "emptyChunkVector.hpp"
#include "pointerHandle.hpp"

/**
 * @brief Global namespace for CoDiPack - Code Differentiation Package
 */
namespace codi {

  /**
   * @brief A vector which manages chunks of data for the taping process.
   *
   * The vector stores an array of data chunks which have all the same size.
   * The data in the chunk can be accessed in a stack like fashion. The user
   * has to check first if enough data is available. The chunk vector will
   * make sure that the current loaded chunk has enough data. The user can then
   * push as many data items as he has reserved on the chunk vector.
   *
   * The read access to the data is provided by the function forEachReverse, which will
   * call the provided function handle on every data item.
   *
   * As some tapes need multiple chunk vectors, the design of the chunk vector reflects
   * this need. The user never knows when a chunk vector pushes a new chunk on the stack
   * but the users needs this information in order to know which range he can evaluate.
   * Therefor the chunk vector needs to store the position of the nested chunk vector
   * every time it pushes a new chunk. The second template argument is used for this
   * kind of usage. It gives the chunk vector the means to access the information it
   * needs.
   *
   * @tparam    ChunkData   The data the chunk vector will store.
   * @tparam NestedVector   A nested chunk vector used for position information
   *                          every time a chunk is pushed.
   */
  template<typename ChunkData, typename NestedVector = EmptyChunkVector>
  class ChunkVector {
    public:

      /**
       * @brief Position of the nested vector
       */
      typedef typename NestedVector::Position NestedPosition;

      /**
       * @brief Typedef of the ChunkData for other classes
       */
      typedef ChunkData ChunkType;

      /**
       * @brief Typedef of NestedVector for other classes
       */
      typedef NestedVector NestedVectorType;

      /**
       * @brief Position of this chunk vector.
       *
       * The position also includes the position of the nested vector,
       * such that the full position of all the chunk vectors
       * is available to the user.
       */
      struct Position {
        size_t chunk; /**< Index of the chunk */
        size_t data;  /**< Data position in the chunk */

        NestedPosition inner; /**< Position of the nested chunk vector */

        /**
         * @brief Default constructor is needed if this position is used as an inner position.
         */
        Position() :
          chunk(0),
          data(0),
          inner() {}

        /**
         * @brief Create the full position for all the nested vectors.
         * @param chunk   Index of the current chunk.
         * @param  data   Index of the data in the current chunk.
         * @param inner   Position of the nested vector.
         */
        Position(const size_t& chunk, const size_t& data, const NestedPosition& inner) :
          chunk(chunk),
          data(data),
          inner(inner) {}

        /**
         * @brief Compares first the inner position and then the own data.
         * @param[in] o  The other position.
         * @return False if the inner position and the own data are not equal.
         */
        bool operator != (const Position& o) {
          return this->inner != o.inner || chunk != o.chunk || data != o.data;
        }

        /**
         * @brief Compares first the inner position and then the own data.
         * @param[in] o  The other position.
         * @return True if the inner position and the own data are not equal.
         */
        bool operator == (const Position& o) {
          return this->inner == o.inner && chunk == o.chunk && data == o.data;
        }
    };

  private:
    std::vector<ChunkData* > chunks; /**< Array of the chunks */

    /**
     * @brief Array of nested positions for the chunks.
     *
     * This vector contains the position of the nested chunk vector. The
     * position is updated if the corresponding chunk in the chunks array
     * is loaded.
     */
    std::vector<NestedPosition> positions;

    ChunkData* curChunk; /**< The current loaded chunk. This will never be NULL */
    size_t curChunkIndex; /**< Index of the chunk which loaded. */

    size_t chunkSize; /**< Global size of the chunks. If this size is set all the chunks are set to this size. */

    NestedVector* nested; /**< Pointer to the nested vector. */

  public:

    /**
     * @brief Creates one chunk and loads it.
     * @param chunkSize   The size for the chunks.
     * @param    nested   The nested chunk vector.
     */
    ChunkVector(const size_t& chunkSize, NestedVector* nested) :
      chunks(),
      positions(),
      curChunk(NULL),
      curChunkIndex(0),
      chunkSize(chunkSize),
      nested(NULL)
    {
      setNested(nested);
    }

    /**
     * @brief Initializes the data structures without touching the nested vector.
     *
     * The method setNested needs to be called to finalize the initialization.
     *
     * @param chunkSize   The size for the chunks.
     */
    ChunkVector(const size_t& chunkSize) :
      chunks(),
      positions(),
      curChunk(NULL),
      curChunkIndex(0),
      chunkSize(chunkSize),
      nested(NULL)
    {}

    /**
     * @brief Deletes all chunks.
     */
    ~ChunkVector() {
      for(size_t i = 0; i < chunks.size(); ++i) {
        delete chunks[i];
      }
    }

    /**
     * @brief Initialize the nested vector.
     *
     * Can only be called once.
     *
     * @param[in,out] v  The nested vector.
     */
    void setNested(NestedVector* v) {
      // Set nested is only called once during the initialization.
      codiAssert(this->nested == NULL);

      this->nested = v;

      curChunk = new ChunkData(chunkSize);
      chunks.push_back(curChunk);
      positions.push_back(nested->getZeroPosition());
    }

    /**
     * @brief Get the nested vector.
     * @return The nested vector.
     */
    NestedVector& getNested() {
      return *nested;
    }

    /**
     * @brief Swap the contents of this chunk vector with the contents of the other
     *        chunk vector.
     *
     * On standard containers the default std::swap method is used.
     * The method is called also on the nested vector.
     *
     * @param[in,out] other  The other chunk vector.
     */
    void swap(ChunkVector<ChunkData, NestedVector>& other) {
      std::swap(chunks, other.chunks);
      std::swap(positions, other.positions);
      std::swap(curChunkIndex, other.curChunkIndex);
      std::swap(chunkSize, other.chunkSize);

      curChunk = chunks[curChunkIndex];
      other.curChunk = other.chunks[other.curChunkIndex];

      nested->swap(*other.nested);

    }

    /**
     * @brief Sets the global chunk size and sets the size of all chunks.
     * @param chunkSize   The new chunk size.
     */
    void setChunkSize(const size_t& chunkSize) {
      this->chunkSize = chunkSize;

      for(size_t i = 0; i < chunks.size(); ++i) {
        chunks[i]->resize(this->chunkSize);
      }
    }

    /**
     * @brief Ensures that enough chunks are allocated so that totalSize data items can be stored.
     * @param totalSize   The number of data items which should be available.
     */
    void resize(const size_t& totalSize) {
      size_t noOfChunks = totalSize / chunkSize;
      if(0 != totalSize % chunkSize) {
        noOfChunks += 1;
      }

      for(size_t i = chunks.size(); i < noOfChunks; ++i) {
        chunks.push_back(new ChunkData(chunkSize));
        positions.push_back(nested->getPosition());
      }
    }

    /**
     * @brief Loads the next chunk.
     *
     * If the current chunk was the last chunk in the array a new chunk is created.
     * Otherwise the old chunk is reset and loaded as current chunk.
     *
     * Always the position of the nested chunk vector is stored.
     */
    CODI_NO_INLINE void nextChunk() {
      curChunk->store();

      curChunkIndex += 1;
      if(chunks.size() == curChunkIndex) {
        curChunk = new ChunkData(chunkSize);
        chunks.push_back(curChunk);
        positions.push_back(nested->getPosition());
      } else {
        curChunk = chunks[curChunkIndex];
        curChunk->reset();
        positions[curChunkIndex] = nested->getPosition();
      }
    }

    /**
     * @brief Resets the chunk vector to the given position.
     *
     * This method will call reset on all chunks which are behind the
     * given position.
     *
     * It calls also reset on the nested chunk vector.
     *
     * @param pos   The position to reset to.
     */
    void reset(const Position& pos) {
      codiAssert(pos.chunk < chunks.size());
      codiAssert(pos.data <= chunkSize);

      for(size_t i = curChunkIndex; i > pos.chunk; --i) {
        chunks[i]->reset();
      }

      curChunk = chunks[pos.chunk];
      curChunk->load();
      curChunk->setUsedSize(pos.data);
      curChunkIndex = pos.chunk;

      nested->reset(pos.inner);
    }

    /**
     * @brief Resets the complete chunk vector.
     */
    void reset() {
      reset(getZeroPosition());
    }

    /**
     * @brief Release all the memory, that the chunk vector has acquired.
     *
     * Reverts the chunk vector into its initial state after the construction.
     * Only the memory of one chunk stays allocated.
     */
    void resetHard() {
      for(size_t i = 1; i < chunks.size(); ++i) {
        delete chunks[i];
      }

      chunks.resize(1);
      positions.resize(1);
      curChunk = chunks[0];
      curChunk->load();

      curChunk->setUsedSize(0);
      curChunkIndex = 0;

      nested->resetHard();
    }

    /**
     * @brief Checks if the current chunk has enough items left.
     *
     * If the chunk has not enough items left, the next chunk is
     * loaded.
     *
     * @param items   The maximum number of items to store.
     */
    CODI_INLINE void reserveItems(const size_t items) {
      codiAssert(items <= chunkSize);

      if(chunkSize < curChunk->getUsedSize() + items) {
        nextChunk();
      }
    }

    /**
     * @brief Sets the data and increases the used chunk data by one.
     *
     * This method should only be called if 'reserveItems' was called
     * beforehand with enough items to accommodate to all calls to this
     * method.
     *
     * @param data  The data set to the current position in the chunk.
     *
     * @tparam Data  The data types for the data to be set.
     */
    template<typename ... Data>
    CODI_INLINE void setDataAndMove(const Data& ... data) {
      // this method should only be called if reserveItems has been called
      curChunk->setDataAndMove(data...);
    }

    /**
     * @brief Sets the provided data pointers to the internal pointers of the
     *        current data.
     *
     * This method should only be called if 'reserveItems' was called
     * beforehand with enough items to accommodate all data that will
     * be written to the pointers.
     *
     * @param[out] pointers  The pointers which are set to the internal data pointers.
     *
     * @tparam Pointers  The data types for the pointers to be set.
     */
    template< typename ... Pointers>
    CODI_INLINE void getDataPointer(Pointers* &... pointers) {
      curChunk->dataPointer(curChunk->usedSize, pointers...);
    }

    /**
     * @brief Advance the count for the data size.
     *
     * This method should only be called if 'reserveItems' was called
     * beforehand with enough items so that this size is smaller or
     * equal.
     *
     * See also #getDataPointer.
     *
     * @param[in] count  The number of items written to the vector.
     */
    CODI_INLINE void addDataSize(size_t count) {
      curChunk->usedSize += count;
    }

    /**
     * @brief The position inside the data of the current chunk.
     * @return The current position in the current chunk.
     */
    CODI_INLINE size_t getChunkPosition() const {
      return curChunk->getUsedSize();
    }

    /**
     * @brief Get the position of the chunk vector and the nested vectors.
     * @return The position of the chunk vector.
     */
    CODI_INLINE Position getPosition() const {
      return Position(curChunkIndex, curChunk->getUsedSize(), nested->getPosition());
    }

    /**
     * @brief Get the zero position of the chunk vector and the nested vectors.
     * @return The zero position of the chunk vector.
     */
    CODI_INLINE Position getZeroPosition() const {
      return Position(0, 0, nested->getZeroPosition());
    }

    /**
     * @brief Get the number of currently allocated chunks.
     * @return The number of currently allocated chunks.
     */
    CODI_INLINE int getNumChunks() const {
      return chunks.size();
    }

    /**
     * @brief Get the chunk size.
     * @return The chunk size.
     */
    CODI_INLINE size_t getChunkSize() const {
      return chunkSize;
    }

    /**
     * @brief Get the total number of data items used.
     * @return The number of data items used in all chunks.
     */
    CODI_INLINE size_t getDataSize() const {
      size_t size = 0;
      for(size_t i = 0; i < chunks.size(); ++i) {
        size += chunks[i]->getUsedSize();
      }

      return size;
    }

  private:
    /**
     * @brief Iterates over the data entries in the chunk.
     *
     * Iterates of the data entries and calls the function object with each data item.
     *
     * It has to hold start >= end.
     *
     * @param chunkPos  The position of the chunk.
     * @param    start  The starting point inside the data of the chunk.
     * @param      end  The end point inside the data of the chunk.
     * @param function  The function called for each data entry.
     * @param     args  Additional arguments for the function.
     *
     * @tparam  Args  The data types for the arguments.
     */
    template<typename FunctionObject, typename ... Args>
    CODI_INLINE void forEachDataReverse(const size_t& chunkPos, const size_t& start, const size_t& end, FunctionObject& function, Args&&... args) {
      codiAssert(start >= end);
      codiAssert(chunkPos < chunks.size());

      PointerHandle<ChunkType> pHandle;

      // we do not initialize dataPos with start - 1 since the type can be unsigned
      for(size_t dataPos = start; dataPos > end; /* decrement is done inside the loop */) {
        --dataPos; // decrement of loop variable

        pHandle.setPointers(dataPos, chunks[chunkPos]);
        pHandle.call(function, std::forward<Args>(args)...);
      }
    }

    /**
     * @brief Iterates over the data entries in the chunk.
     *
     * Iterates of the data entries and calls the function object with each data item.
     *
     * It has to hold start <= end.
     *
     * @param chunkPos  The position of the chunk.
     * @param    start  The starting point inside the data of the chunk.
     * @param      end  The end point inside the data of the chunk.
     * @param function  The function called for each data entry.
     * @param     args  Additional arguments for the function.
     *
     * @tparam  Args  The data types for the arguments.
     */
    template<typename FunctionObject, typename ... Args>
    CODI_INLINE void forEachDataForward(const size_t& chunkPos, const size_t& start, const size_t& end, FunctionObject& function, Args&&... args) {
      codiAssert(start <= end);
      codiAssert(chunkPos < chunks.size());

      PointerHandle<ChunkType> pHandle;

      for(size_t dataPos = start; dataPos < end; dataPos += 1) {
        pHandle.setPointers(dataPos, chunks[chunkPos]);
        pHandle.call(function, std::forward<Args>(args)...);
      }
    }

  public:

    /**
     * @brief Iterates over all data entries in the given range
     *
     * Iterates of the data entries and calls the function object with each data item.
     *
     * It has to hold start >= end.
     *
     * @param    start  The starting point of the range.
     * @param      end  The end point of the range.
     * @param function  The function called for each data entry.
     * @param     args  Additional arguments for the function.
     *
     * @tparam  Args  The data types for the arguments.
     */
    template<typename FunctionObject, typename ... Args>
    CODI_INLINE void forEachReverse(const Position& start, const Position& end, FunctionObject& function, Args&&... args) {
      codiAssert(start.chunk > end.chunk || (start.chunk == end.chunk && start.data >= end.data));
      codiAssert(start.chunk < chunks.size());

      size_t dataStart = start.data;
      for(size_t chunkPos = start.chunk; chunkPos > end.chunk; /* decrement is done inside the loop */) {

        forEachDataReverse(chunkPos, dataStart, 0, function, std::forward<Args>(args)...);

        dataStart = chunks[--chunkPos]->getUsedSize(); // decrement of loop variable

      }

      forEachDataReverse(end.chunk, dataStart, end.data, function, std::forward<Args>(args)...);
    }

    /**
     * @brief Iterates over all data entries in the given range in a forward loop.
     *
     * Iterates of the data entries and calls the function object with each data item.
     *
     * It has to hold start <= end.
     *
     * @param    start  The starting point of the range.
     * @param      end  The end point of the range.
     * @param function  The function called for each data entry.
     * @param     args  Additional arguments for the function.
     *
     * @tparam  Args  The data types for the arguments.
     */
    template<typename FunctionObject, typename ... Args>
    CODI_INLINE void forEachForward(const Position& start, const Position& end, FunctionObject& function, Args&&... args) {
      codiAssert(start.chunk < end.chunk || (start.chunk == end.chunk && start.data <= end.data));
      codiAssert(end.chunk < chunks.size());

      size_t dataStart = start.data;
      for(size_t chunkPos = start.chunk; chunkPos < end.chunk; chunkPos += 1) {

        forEachDataForward(chunkPos, dataStart, chunks[chunkPos]->getUsedSize(), function, std::forward<Args>(args)...);

        dataStart = 0;

      }

      forEachDataForward(end.chunk, dataStart, end.data, function, std::forward<Args>(args)...);
    }

    /**
     * @brief Iterates over all chunks.
     *
     * Iterates over all chunks of the vector. If the recursive argument is given,
     * the iteration continuous with the chunks from the nested vector.
     *
     * The function object will be called with the chunk as the first argument, followed by the given arguments args.
     *
     * @param  function  The function called for each chunk.
     * @param recursive  If also the chunks of the nested vectors should be iterated.
     * @param      args  The pointers are used as the arguments for the function.
     *
     * @tparam  Args  The data types for the arguments of the function.
     */
    template<typename FunctionObject, typename ... Args>
    CODI_INLINE void forEachChunkForward(FunctionObject& function, bool recursive, Args&&... args) {

      for(size_t chunkPos = 0; chunkPos < chunks.size(); ++chunkPos) {

        function(chunks[chunkPos], std::forward<Args>(args)...);
      }

      if(recursive) {
        nested->forEachChunkForward(function, recursive, std::forward<Args>(args)...);
      }
    }

    /**
     * @brief Reverse stack evaluation of the tape.
     *
     * All pointers to the data items are created and given with the start and end position for
     * the interpretation range to the next vector. The last vector will call the provided function.
     *
     * The function is called several times for each valid range described by the chunks of the nested
     * vectors. The function has to modify the dataPos given for each chunk vector such that it is reduced
     * to end position for the interpretation.
     *
     * The function call is
     * \code{.cpp}
     * func(start.nested, end.nested, <other arguments>,
     *      startDataPos, endDataPos, pointerChunkItem1, pointerChunkItem2, etc.);
     * \endcode
     *
     * Debug checks will ensure this behaviour.
     *
     * It has to hold start >= end.
     *
     * @param    start  The start point for the stack interpretation.
     * @param      end  The end point for the stack interpretation.
     * @param function  The function called for each valid range.
     * @param     args  Pointers and ranges from other chunks vectors and additional arguments for the
     *                  function.
     *
     * @tparam  Args  The data types for the arguments.
     */
    template<typename Function, typename ... Args>
    CODI_INLINE void evaluateReverse(const Position& start, const Position& end,const Function& function,
                                     Args&&... args) {
      PointerHandle<ChunkType> pHandle;

      size_t dataPos = start.data;
      NestedPosition curInnerPos = start.inner;
      for(size_t curChunk = start.chunk; curChunk > end.chunk; --curChunk) {

        pHandle.setPointers(0, chunks[curChunk]);

        NestedPosition endInnerPos = positions[curChunk];
        pHandle.callNestedReverse(nested, curInnerPos, endInnerPos, function, std::forward<Args>(args)..., dataPos, 0);

        codiAssert(dataPos == 0); // after a full chunk is evaluated, the data position needs to be zero

        curInnerPos = endInnerPos;

        dataPos = chunks[curChunk - 1]->getUsedSize();
      }

      // Iterate over the reminder also covers the case if the start chunk and end chunk are the same
      pHandle.setPointers(0, chunks[end.chunk]);
      pHandle.callNestedReverse(nested, curInnerPos, end.inner, function, std::forward<Args>(args)..., dataPos, end.data);

      codiAssert(dataPos == end.data); // after the last chunk is evaluated, the data position needs to be at the end position
    }

    /**
     * @brief Forward stack evaluation of the tape.
     *
     * All pointers to the data items are created and given with the start and end position for
     * the interpretation range to the next vector. The last vector will call the provided function.
     *
     * The function is called several times for each valid range described by the chunks of the nested
     * vectors. The function has to modify the dataPos given for each chunk vector such that it is increased
     * to end position for the interpretation.
     *
     * The function call is
     * \code{.cpp}
     * func(start.nested, end.nested, <other arguments>,
     *      startDataPos, endDataPos, pointerChunkItem1, pointerChunkItem2, etc.);
     * \endcode
     *
     * Debug checks will ensure this behaviour.
     *
     * It has to hold start >= end.
     *
     * @param    start  The start point for the stack interpretation.
     * @param      end  The end point for the stack interpretation.
     * @param function  The function called for each valid range.
     * @param     args  Pointers and ranges from other chunks vectors and additional arguments for the
     *                  function.
     *
     * @tparam  Args  The data types for the arguments.
     */
    template<typename Function, typename ... Args>
    CODI_INLINE void evaluateForward(const Position& start, const Position& end,const Function& function,
                                     Args&&... args) {
      PointerHandle<ChunkType> pHandle;

      size_t dataPos = start.data;
      NestedPosition curInnerPos = start.inner;
      for(size_t curChunk = start.chunk; curChunk < end.chunk; ++curChunk) {

        pHandle.setPointers(0, chunks[curChunk]);

        NestedPosition endInnerPos = positions[curChunk + 1];
        pHandle.callNestedForward(nested, curInnerPos, endInnerPos, function, std::forward<Args>(args)..., dataPos, chunks[curChunk]->getUsedSize());

        // After a full chunk is evaluated, the data position needs to be at the end of the chunk
        codiAssert(dataPos == chunks[curChunk]->getUsedSize());

        curInnerPos = endInnerPos;

        dataPos = 0;
      }

      // Iterate over the reminder also covers the case if the start chunk and end chunk are the same
      pHandle.setPointers(0, chunks[end.chunk]);
      pHandle.callNestedForward(nested, curInnerPos, end.inner, function, std::forward<Args>(args)..., dataPos, end.data);

      codiAssert(dataPos == end.data); // after the last chunk is evaluated, the data position needs to be at the end position
    }
  };
}
