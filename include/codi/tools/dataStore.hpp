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
#include <algorithm>
#include <vector>

/**
 * @brief Global namespace for CoDiPack - Code Differentiation Package
 */
namespace codi {

  /**
   * @brief Base class the data in the DataStore.
   *
   * It just defines a void pointer to the data. The logic is
   * implemented in the template versions of the class.
   */
  class DataHandleBase {
  public:
    /** @brief The pointer to the data. */
    void* data;

    /** @brief Virtual delete function. */
    virtual ~DataHandleBase() {}

    /**
     * @brief The clone method is used to make a deep copy of the DataStore.
     * @return A new data object with the same data.
     */
    virtual DataHandleBase* clone() = 0;
  };

  /**
   * @brief Default template implementation for data in the data store.
   *
   * The implementation casts the void pointer from the base class to
   * the given type.
   *
   * @tparam The type of the data.
   */
  template<typename Type>
  class DataHandle : public DataHandleBase {
  public:
    /**
     * @brief The constructor will copy the given data.
     * @param[in] value The data for this object.
     */
    explicit DataHandle(const Type& value) {
      data = (void*) new Type(value);
    }

    /**
     * @brief Frees the data object.
     */
    ~DataHandle() {
      Type* pointer = (Type*) data;
      delete pointer;
    }

    /** @brief Constructs a new data handle with the same data */
    DataHandleBase* clone() { return new DataHandle<Type>(*((Type*) data)); }
  };

  /**
   * @brief Default template implementation for array data in the data store.
   *
   * The implementation casts the void pointer from the base class to
   * the given type.
   *
   * @tparam The type of the data in the array.
   */
  template<typename Type>
  class DataHandleArray : public DataHandleBase {
  private:
    /** @brief The size of the array. */
    int size;
  public:
    /**
     * @brief The constructor will copy the given data.
     * @param[in] value The data for this object.
     * @param[in]  size The size of the data array.
     */
    DataHandleArray(const Type* value, int size) {
      data = (void*) new Type[size];
      this->size = size;
      std::copy(value, &value[size], (Type*) data);
    }

    /**
     * @brief Frees the data object.
     */
    ~DataHandleArray() {
      Type* pointer = (Type*) data;
      delete[] pointer;
    }

    /** @brief Constructs a new data handle with the same data */
    DataHandleBase* clone() { return new DataHandleArray<Type>(*((Type*) data), size); }
  };

  /**
   * @brief The class is used to store data for external functions.
   *
   * The DataStore saves the data which is given to it. All the data is copied and
   * no pointers to the original data are stored.
   *
   * The store can also store pointers to the data. Then the users has to take
   * care that the data is not modified and valid in the instance were the data is
   * used.
   *
   * The order in which the data is read and written to the data store has to be the same.
   * Otherwise the behaviour is not guaranteed. When all data has been read from the data store
   * it will reset such that the next call to a getData function will return the first item.
   */
  class DataStore {
  private:

    /** @brief The vector with the data handles. */
    std::vector<DataHandleBase*> store;

    /** @brief The current position in the data vector.*/
    size_t storePos;

  public:

    /**
     * @brief Create an empty data store.
     */
    DataStore() :
      storePos(0) { }

    /**
     * @brief Deletes all the data handles.
     */
    ~DataStore() {
      clear();
    }

    /**
     * @brief Deletes all the data handles.
     */
    void clear() {
      for(size_t i = 0; i < store.size(); ++i) {
        delete store[i];
      }
      store.clear();
    }

    /**
     * @brief Copy constructor. Creates a deep copy of the data in the data store.
     * @param[in] other The data store which is cloned.
     */
    DataStore(const DataStore& other) {
      for(size_t i = 0; i < other.store.size(); ++i) {
        store.push_back(other.store[i]->clone());
      }
      storePos = other.storePos;
    }

    /**
     * @brief Copy operator. Creates a deep copy of the data in the data store.
     * @param[in] other The data store which is cloned.
     */
    DataStore& operator=(const DataStore& other) {
      this->clear();
      for(size_t i = 0; i < other.store.size(); ++i) {
        store.push_back(other.store[i]->clone());
      }
      storePos = other.storePos;

      return *this;
    }

    /**
     * @brief Add data to the data store.
     *
     * @param[in] value  The data which is stored.
     * @return The position of the added data.
     *
     * @tparam Type The type of the stored data.
     */
    template<typename Type>
    size_t addData(const Type& value) {
      store.push_back(new DataHandle<Type>(value));
      return store.size() - 1;
    }

    /**
     * @brief Add array data to the data store.
     *
     * @param[in] value  The data pointer to the array.
     * @param[in]  size  The size of the stored data.
     * @return The position of the added data.
     *
     * @tparam Type The type of the stored array data.
     */
    template<typename Type>
    size_t addData(const Type* value, const int size) {
      store.push_back(new DataHandleArray<Type>(value, size));
      return store.size() - 1;
    }

    /**
     * @brief Get data from the data store.
     *
     * The call to this function has to correspond to
     * the write call which stored the actual data.
     *
     * @param[out] value  The data will be set to this value.
     *
     * @tparam Type The type of the extracted data.
     */
    template<typename Type>
    void getData(Type& value) {
      getData<Type>(&value, 1);
    }

    /**
     * @brief Get data from the data store, without copying it.
     *
     * The call to this function has to correspond to
     * the write call which stored the actual data.
     *
     * @return The const reference to the data stored.
     *
     * @tparam Type The type of the extracted data.
     */
    template<typename Type>
    const Type& getData() {
      Type* data = nextStore<Type>();
      return *data;
    }

    /**
     * @brief Get data from the data store.
     *
     * The call to this function has to correspond to
     * the write call which stored the actual data.
     *
     * @param[out] value  The data will be copied to the value.
     * @param[in]   size  The size of the data array.
     *
     * @tparam Type The type of the extracted data.
     */
    template<typename Type>
    void getData(Type* value, const int size) {
      Type* convPointer = nextStore<Type>();

      std::copy(convPointer, &convPointer[size], value);
    }

    /**
     * @brief Get array data from the data store, without copying it.
     *
     * The call to this function has to correspond to
     * the write call which stored the actual data.
     *
     * @return The const pointer to the data stored.
     *
     * @tparam Type The type of the extracted data.
     */
    template<typename Type>
    const Type* getDataArray() {
      Type* data = nextStore<Type>();
      return *data;
    }

    /**
     * @brief Get data from the data store with the index from the add function.
     *
     * @param[out] value  The data will be set to this value.
     * @param[in]    pos  The position for the data. This needs to be the index returned by the add function.
     *
     * @tparam Type The type of the extracted data.
     */
    template<typename Type>
    void getDataByIndex(Type& value, size_t pos) {
      getDataArrayByIndex<Type>(&value, 1, pos);
    }

    /**
     * @brief Get data from the data store with the index from the add function.
     *
     * @param[out] value  The data will be copied to the value.
     * @param[in]   size  The size of the data array.
     * @param[in]    pos  The position for the data. This needs to be the index returned by the add function.
     *
     * @tparam Type The type of the extracted data.
     */
    template<typename Type>
    void getDataArrayByIndex(Type* value, const int size, size_t pos) {
      Type* convPointer = getStore<Type>(pos);

      std::copy(convPointer, &convPointer[size], value);
    }

    /**
     * @brief Restart the reading process.
     */
    void resetPos() {
      storePos = 0;
    }

  private:

    /**
     * Gets a specific data item from the data vector.
     *
     * @param[in] pos  The position of the item.
     * @return The data item for the specified position.
     *
     * @tparam Type The type of the extracted data.
     */
    template<typename Type>
    Type* getStore(size_t pos) {
      return (Type*) store[pos]->data;
    }

    /**
     * @brief Gets the next data item out of the data vector.
     *
     * The call to this function has to correspond to
     * the write call which stored the actual data.
     *
     * @return The pointer to the data stored.
     *
     * @tparam Type The type of the extracted data.
     */
    template<typename Type>
    Type* nextStore() {
      Type* pointer = getStore<Type>(storePos);
      storePos += 1;
      if(storePos >= store.size()) {
        storePos = 0;
      }

      return pointer;
    }
  };
}
