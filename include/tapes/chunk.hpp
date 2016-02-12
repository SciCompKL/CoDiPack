/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015 Chair for Scientific Computing (SciComp), TU Kaiserslautern
 * Homepage: http://www.scicomp.uni-kl.de
 * Contact:  Prof. Nicolas R. Gauger (codi@scicomp.uni-kl.de)
 *
 * Lead developers: Max Sagebaum, Tim Albring (SciComp, TU Kaiserslautern)
 *
 * This file is part of CoDiPack (http://www.scicomp.uni-kl.de/software/codi).
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
 * Authors: Max Sagebaum, Tim Albring, (SciComp, TU Kaiserslautern)
 */
#pragma once

#include <cstddef>
#include <new>
#include <string.h>

#include "../configure.h"

/**
 * @brief Global namespace for CoDiPack - Code Differentiation Package
 */
namespace codi {

  /**
   * @brief The basic interface for data chunks.
   *
   * The chunk interface is the basic interface all chunks have to implement.
   * The interface already provides the basic facilities for the
   * size of the data inside the chunk.
   */
  struct ChunkInterface {

    size_t size; /**< Size of the allocated data */
    size_t usedSize; /**< Number of used items in the data array */

    /**
     * @brief Create a chunk with the given size.
     * @param size  The size of the data in the chunk.
     */
    ChunkInterface(const size_t& size) :
      size(size),
      usedSize(0) {}

    /**
     * @brief Get the maximum size of the chunk
     * @return The maximum number of items.
     */
    inline size_t getSize() const {
      return size;
    }

    /**
     * @brief Get the number of used items.
     * @return The number of used items.
     */
    inline size_t getUsedSize() const {
      return usedSize;
    }

    /**
     * @brief Get the number of free items.
     * @return The number of free items.
     */
    inline size_t getUnusedSize() const {
      return size - usedSize;
    }

    /**
     * @brief Fully reset the data in this chunk.
     */
    inline void reset() {
      usedSize = 0;
    }

    /**
     * @brief Set the number of used items in this chunk.
     * @param usage   The number of used items.
     */
    inline void setUsedSize(const size_t& usage) {
      usedSize = usage;
    }

    /**
     * @brief Store the data of the chunk.
     *
     * This method is called when the data of a chunk is no
     * longer directly needed and can be stored somewhere else.
     */
    void store() {}

    /**
     * @brief Load the data of the chunk.
     *
     * This method is called when the data of a chunk is needed by the
     * evaluation process.
     */
    void load() {}
  };

  /**
   * @brief Chunk with one data array.
   *
   * This chunk contains one data array which is stored in memory.
   *
   * @tparam Data   The type of the stored data. This type has to be a POD type as we use new, free and memset on it.
   */
  template<typename Data>
  struct Chunk1 : public ChunkInterface {
    Data* data; /**< The data of the chunk */

    /**
     * @brief Creates the data of the chunk.
     *
     * The data of the chunk is set with memset to zero.
     *
     * @param size The size of the data in the chunk.
     */
    Chunk1(const size_t& size) : ChunkInterface(size) {
      data = (Data*)malloc(sizeof(Data) * size);
      if(UseMemsetInChunks) {
        memset(data, 0, sizeof(Data) * size);
      }
    }

    /**
     * @brief Deletes the data array
     */
    ~Chunk1() {
      free(data);
      data = NULL;
    }

    /**
     * @brief Set the size of the array.
     * @param size  The new size of the array.
     */
    void resize(const size_t &size) {
      this->~Chunk1();
      new (this) Chunk1(size);
    }

    /**
     * @brief Set the data values to the current position and increment the used size.
     * @param values  The values which are set to the data.
     */
    inline void setDataAndMove(const Data& values) {
      assert(getUnusedSize() != 0);
      data[usedSize] = values;
      ++usedSize;
    }

    /**
     * @brief Returns a pointer to the data array at the given position.
     * @param index       The index in the data array.
     * @param dataPointer Pointer that is set to the internal data pointer.
     */
    inline void dataPointer(const size_t& index, Data* &dataPointer) {
      assert(index <= ChunkInterface::size);

      dataPointer = &data[index];
    }
  };

  /**
   * @brief Chunk with two data arrays.
   *
   * This chunk contains two data arrays which are stored in memory.
   *
   * @tparam Data1   The first type of the stored data. This type has to be a POD type as we use new, free and memset on it.
   * @tparam Data2   The second type of the stored data. This type has to be a POD type as we use new, free and memset on it.
   */
  template<typename Data1, typename Data2>
  struct Chunk2 : public ChunkInterface {
    Data1* data1; /**< First data item of the chunk */
    Data2* data2; /**< Second data item of the chunk */

    /**
     * @brief Creates the data of the chunk.
     *
     * The data of the chunk is set with memset to zero.
     *
     * @param size The size of the data in the chunk.
     */
    Chunk2(const size_t& size) : ChunkInterface(size) {
      data1 = (Data1*)malloc(sizeof(Data1) * size);
      data2 = (Data2*)malloc(sizeof(Data2) * size);
      if(UseMemsetInChunks) {
        memset(data1, 0, sizeof(Data1) * size);
        memset(data2, 0, sizeof(Data2) * size);
      }
    }

    /**
     * @brief Deletes the data arrays
     */
    ~Chunk2() {
      free(data1);
      free(data2);
      data1 = NULL;
      data2 = NULL;
    }

    /**
     * @brief Set the size of the arrays.
     * @param size  The new size of the arrays.
     */
    void resize(const size_t &size) {
      this->~Chunk2();
      new (this) Chunk2(size);
    }

    /**
     * @brief Set the data values to the current position and increment the used size.
     * @param value1  The value for the first data array.
     * @param value2  The value for the second data array.
     */
    inline void setDataAndMove(const Data1& value1, const Data2& value2) {
      assert(getUnusedSize() != 0);
      data1[usedSize] = value1;
      data2[usedSize] = value2;
      ++usedSize;
    }

    /**
     * @brief Returns a pointer to the data array at the given position.
     * @param index   The index in the data array.
     * @return A pointer to the data.
     */
    inline void dataPointer(const size_t& index, Data1* &pointer1, Data2* &pointer2) {
      assert(index <= ChunkInterface::size);
      pointer1 = &data1[index];
      pointer2 = &data2[index];
    }
  };

}
