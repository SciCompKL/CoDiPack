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
   * @brief A vector that creates only one chunk and presents it like a ChunkVector.
   *
   * The vector stores one data chunk.
   * The data in the chunk can be accessed in a stack like fashion. The user
   * has to ensure that enough data is present. All the usual checks with reserveItems
   * are only performed in codiAssert statements.
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
   */
  template<typename ChunkData, typename NestedVector = EmptyChunkVector>
  class SingleChunkVector {
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
     * @brief Position of the nested vector
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
      size_t data;  /**< Data position in the chunk (Will always be zero. Required for compatibility to ChunkVector. */

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
       * @param  data   Index of the data in the current chunk.
       * @param inner   Position of the nested vector.
       */
      Position(const size_t& data, const NestedPosition& inner) :
        chunk(0),
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

    ChunkData chunk; /**< The data chunk. */

    NestedVector* nested; /**< Pointer to the nested vector. */

  public:

    /**
     * @brief Creates one chunk and loads it.
     * @param chunkSize   The size for the chunks.
     * @param    nested   The nested chunk vector.
     */
    SingleChunkVector(const size_t& chunkSize, NestedVector* nested) :
      chunk(chunkSize),
      nested(nested)
    {}

    /**
     * @brief Initializes the data structures without touching the nested vector.
     *
     * The method setNested needs to be called to finalize the initialization.
     *
     * @param chunkSize   The size for the chunks.
     */
    SingleChunkVector(const size_t& chunkSize) :
      chunk(chunkSize),
      nested(NULL)
    {}

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
    void swap(SingleChunkVector<ChunkData, NestedVector>& other) {
      chunk.swap(other.chunk);

      nested->swap(*other.nested);

    }

    /**
     * @brief Sets the size of the chunk.
     * @param chunkSize   The new chunk size.
     */
    void setChunkSize(const size_t& chunkSize) {
       chunk.resize(chunkSize);
    }

    /**
     * @brief Ensures that the chunk has enough size available.
     * @param totalSize   The number of data items which should be available.
     */
    void resize(const size_t& totalSize) {
      setChunkSize(totalSize);
    }

    /**
     * @brief Resets the chunk vector to the given position.
     *
     * This method will set the used size of the chunk to the one in the position.
     *
     * It also calls reset on the nested chunk vector.
     *
     * @param pos   The position to reset to.
     */
    void reset(const Position& pos) {
      codiAssert(pos.data < chunk.getSize());

      chunk.setUsedSize(pos.data);

      nested->reset(pos.inner);
    }

    /**
     * @brief Resets the complete chunk vector.
     */
    void reset() {
      reset(getZeroPosition());
    }

    /**
     * @brief Only calls the nested vector.
     */
    void resetHard() {
      nested->resetHard();
    }

    /**
     * @brief Performs no check only does an codiAssert.
     *
     * @param items   The maximum number of items to store.
     */
    CODI_INLINE void reserveItems(const size_t items) const {
      codiAssert(chunk.getUsedSize() + items < chunk.getSize());
    }

    /**
     * @brief Sets the data and increases the used chunk data by one.
     *
     * This method should only be called if 'reserveItems' was called
     * beforehand with enough items to accommodate to all calls to this
     * method.
     *
     * In this implementation the user has to ensure that there is enough
     * space allocated up front.
     *
     * @param data  The data set to the current position in the chunk.
     *
     * @tparam Data  The data types for the data to be set.
     */
    template<typename ... Data>
    CODI_INLINE void setDataAndMove(const Data& ... data) {
      // this method should only be called if reserveItems has been called
      chunk.setDataAndMove(data...);
    }

    /**
     * @brief The position inside the data of the current chunk.
     * @return The current position in the current chunk.
     */
    CODI_INLINE size_t getChunkPosition() const {
      return chunk.getUsedSize();
    }

    /**
     * @brief Get the position of the chunk vector and the nested vectors.
     * @return The position of the chunk vector.
     */
    CODI_INLINE Position getPosition() const {
      return Position(chunk.getUsedSize(), nested->getPosition());
    }

    /**
     * @brief Get the zero position of the chunk vector and the nested vectors.
     * @return The zero position of the chunk vector.
     */
    CODI_INLINE Position getZeroPosition() const {
      return Position(0, nested->getZeroPosition());
    }

    /**
     * @brief Get the number of currently allocated chunks.
     * @return The number of currently allocated chunks.
     */
    CODI_INLINE int getNumChunks() const {
      return 1;
    }

    /**
     * @brief Get the chunk size.
     * @return The chunk size.
     */
    CODI_INLINE size_t getChunkSize() const {
      return chunk.getSize();
    }

    /**
     * @brief Get the total number of data items used.
     * @return The number of data items used in all chunks.
     */
    CODI_INLINE size_t getDataSize() const {
      return chunk.getUsedSize();
    }

  private:
    /**
     * @brief Iterates over the data entries in the chunk.
     *
     * Iterates of the data entries and calls the function object with each data item.
     *
     * It has to hold start >= end.
     *
     * @param    start  The starting point inside the data of the chunk.
     * @param      end  The end point inside the data of the chunk.
     * @param function  The function called for each data entry.
     * @param     args  Additional arguments for the function.
     *
     * @tparam  Args  The data types for the arguments.
     */
    template<typename FunctionObject, typename ... Args>
    CODI_INLINE void forEachDataReverse(const size_t& start, const size_t& end, FunctionObject& function, Args&&... args) {
      codiAssert(start >= end);

      PointerHandle<ChunkType> pHandle;

      // we do not initialize dataPos with start - 1 because the type can be unsigned
      for(size_t dataPos = start; dataPos > end; /* decrement is done inside the loop */) {
        --dataPos; // decrement of loop variable

        pHandle.setPointers(dataPos, &chunk);
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
     * @param    start  The starting point inside the data of the chunk.
     * @param      end  The end point inside the data of the chunk.
     * @param function  The function called for each data entry.
     * @param     args  Additional arguments for the function.
     *
     * @tparam  Args  The data types for the arguments.
     */
    template<typename FunctionObject, typename ... Args>
    CODI_INLINE void forEachDataForward(const size_t& start, const size_t& end, FunctionObject& function, Args&&... args) {
      codiAssert(start <= end);

      PointerHandle<ChunkType> pHandle;

      // we do not initialize dataPos with start - 1 because the type can be unsigned
      for(size_t dataPos = start; dataPos < end; dataPos += 1) {
        pHandle.setPointers(dataPos, &chunk);
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
     * @param     args  Additional arguments for the function
     *
     * @tparam  Args  The data types for the arguments.
     */
    template<typename FunctionObject, typename ... Args>
    CODI_INLINE void forEachReverse(const Position& start, const Position& end, FunctionObject& function, Args &&... args) {
      codiAssert(start.chunk == 0);
      codiAssert(end.chunk == 0);
      codiAssert(start.data >= end.data);
      codiAssert(start.data <= chunk.getSize());

      forEachDataReverse(start.data, end.data, function, std::forward<Args>(args)...);
    }

    /**
     * @brief Iterates over all data entries in the given range
     *
     * Iterates of the data entries and calls the function object with each data item.
     *
     * It has to hold start <= end.
     *
     * @param    start  The starting point of the range.
     * @param      end  The end point of the range.
     * @param function  The function called for each data entry.
     * @param     args  Additional arguments for the function
     *
     * @tparam  Args  The data types for the arguments.
     */
    template<typename FunctionObject, typename ... Args>
    CODI_INLINE void forEachForward(const Position& start, const Position& end, FunctionObject& function, Args &&... args) {
      codiAssert(start.chunk == 0);
      codiAssert(end.chunk == 0);
      codiAssert(start.data <= end.data);
      codiAssert(end.data <= chunk.getSize());

      forEachDataForward(start.data, end.data, function, std::forward<Args>(args)...);
    }

    /**
     * @brief Iterates over the chunk of the vector.
     *
     * If the recursive argument is given the iteration continuous with the chunks from the nested vector.
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
    CODI_INLINE void forEachChunkForward(FunctionObject& function, bool recursive, Args &&... args) {

      function(&chunk, args...);

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
      pHandle.setPointers(0, &chunk);
      pHandle.callNestedReverse(nested, start.inner, end.inner, function, std::forward<Args>(args)..., dataPos, end.data);
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
      pHandle.setPointers(0, &chunk);
      pHandle.callNestedForward(nested, start.inner, end.inner, function, std::forward<Args>(args)..., dataPos, end.data);
    }
  };
}
