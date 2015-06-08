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
#include <algorithm>
#include <vector>

namespace codi {

  class DeleteHandleBase {
  public:
    void* data;

    virtual ~DeleteHandleBase() { };

    virtual DeleteHandleBase* clone() = 0;
  };

  template<typename Type>
  class DeleteHandle : public DeleteHandleBase {
  public:
    DeleteHandle(const Type& value) {
      data = (void*) new Type(value);
    }

    ~DeleteHandle() {
      Type* pointer = (Type*) data;
      delete pointer;
    }

    DeleteHandleBase* clone() { return new DeleteHandle<Type>(*((Type*) data)); }
  };

  template<typename Type>
  class DeleteHandleArray : public DeleteHandleBase {
  private:
    int size;
  public:
    DeleteHandleArray(const Type* value, int size) {
      data = (void*) new Type[size];
      this->size = size;
      std::copy(value, &value[size], (Type*) data);
    }

    ~DeleteHandleArray() {
      Type* pointer = (Type*) data;
      delete[] pointer;
    }

    DeleteHandleBase* clone() { return new DeleteHandleArray<Type>(*((Type*) data), size); }
  };

  class DataStore {
  private:
    std::vector<DeleteHandleBase*> store;

    size_t storePos;

  public:

    DataStore() :
      storePos(0) { }

    ~DataStore() {
      clear();
    }

    void clear() {
      for(size_t i = 0; i < store.size(); ++i) {
        delete store[i];
      }
      store.clear();
    }

    DataStore(const DataStore& other) {
      for(size_t i = 0; i < other.store.size(); ++i) {
        store.push_back(other.store[i]->clone());
      }
      storePos = other.storePos;
    }

    DataStore& operator=(const DataStore& other) {
      this->clear();
      for(size_t i = 0; i < other.store.size(); ++i) {
        store.push_back(other.store[i]->clone());
      }
      storePos = other.storePos;

      return *this;
    }

    template<typename Type>
    void addData(const Type& value) {
      store.push_back(new DeleteHandle<Type>(value));
    }

    template<typename Type>
    void addData(const Type* value, const int size) {
      store.push_back(new DeleteHandleArray<Type>(value, size));
    }

    template<typename Type>
    void getData(Type& value) {
      getData(&value, 1);
    }

    template<typename Type>
    void getData(Type* value, const int size) {
      Type* convPointer = nextStore<Type>();

      std::copy(convPointer, &convPointer[size], value);
    }

    void resetPos() {
      storePos = 0;
    }

  private:
    template<typename Type>
    Type* nextStore() {
      Type* pointer = (Type*) store[storePos++]->data;
      if(storePos >= store.size()) {
        storePos = 0;
      }

      return pointer;
    }
  };
}