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

#include <cstddef>
#include <new>
#include <string.h>
#include <tuple>

namespace codi {

  struct ChunkInterface {

    size_t size;
    size_t usedSize;

    ChunkInterface(const size_t& size) :
      size(size),
      usedSize(0) {}

    ~ChunkInterface() {}

    inline size_t getUsedSize() {
      return usedSize;
    }

    inline size_t getUnusedSize() {
      return size - usedSize;
    }

    inline void reset() {
      usedSize = 0;
    }

    inline size_t increase() {
      return ++usedSize;
    }

    inline void setUsedSize(const size_t& usage) {
      usedSize = usage;
    }

    void store() {}
    void load() {}
  };

  template<typename Data>
  struct Chunk1 : public ChunkInterface {
    typedef std::tuple<Data*> DataPointer;
    typedef std::tuple<Data> DataValues;

    Data* data;

    Chunk1(const size_t& size) : ChunkInterface(size) {
      data = (Data*)malloc(sizeof(Data) * size);
      memset(data, 0, sizeof(Data) * size);
    }

    ~Chunk1() {
      free(data);
      data = NULL;
    }

    void resize(const size_t &size) {
      this->~Chunk1();
      new (this) Chunk1(size);
    }

    inline void setDataAndMove(const DataValues& values) {
      std::tie(data[usedSize]) = values;
      ++usedSize;
    }

    inline std::tuple<Data*> dataPointer(const size_t& index) {
      return std::make_tuple(&data[index]);
    }
  };

  template<typename Data1, typename Data2>
  struct Chunk2 : public ChunkInterface {
    typedef std::tuple<Data1*, Data2*> DataPointer;
    typedef std::tuple<Data1, Data2> DataValues;

    Data1* data1;
    Data2* data2;

    Chunk2(const size_t& size) : ChunkInterface(size) {
      data1 = (Data1*)malloc(sizeof(Data1) * size);
      data2 = (Data2*)malloc(sizeof(Data2) * size);
      memset(data1, 0, sizeof(Data1) * size);
      memset(data2, 0, sizeof(Data2) * size);
    }

    ~Chunk2() {
      free(data1);
      free(data2);
      data1 = NULL;
      data2 = NULL;
    }

    void resize(const size_t &size) {
      this->~Chunk2();
      new (this) Chunk2(size);
    }

    inline void setDataAndMove(const DataValues& values) {
      std::tie(data1[usedSize], data2[usedSize]) = values;
      ++usedSize;
    }

    inline DataPointer dataPointer(const size_t& index) {
      return std::make_tuple(&data1[index], &data2[index]);
    }
  };
}
