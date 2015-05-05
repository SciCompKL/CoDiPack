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

namespace codi {

  struct ChunkInterface {

    size_t size;
    size_t usedSize;

    ChunkInterface(const size_t& size) :
      size(size),
      usedSize(0) {}

    ~ChunkInterface() {}

    virtual void resize(const size_t& size) {
      this->size = size;
    }

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
    Data* data;

    Chunk1(const size_t& size) : ChunkInterface(size) {
      data = new Data[size];
    }

    ~Chunk1() {
      delete [] data;
    }

    virtual void resize(const size_t &size) {
      this->~Chunk1();
      new (this) Chunk1(size);
    }
  };

  template<typename Data1, typename Data2>
  struct Chunk2 : public ChunkInterface {
    Data1* data1;
    Data2* data2;

    Chunk2(const size_t& size) : ChunkInterface(size) {
      data1 = new Data1[size];
      data2 = new Data2[size];
    }

    ~Chunk2() {
      delete [] data1;
      delete [] data2;
    }

    virtual void resize(const size_t &size) {
      this->~Chunk2();
      new (this) Chunk2(size);
    }
  };
}
