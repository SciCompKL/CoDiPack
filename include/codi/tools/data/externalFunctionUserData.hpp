/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2023 Chair for Scientific Computing (SciComp), University of Kaiserslautern-Landau
 * Homepage: http://www.scicomp.uni-kl.de
 * Contact:  Prof. Nicolas R. Gauger (codi@scicomp.uni-kl.de)
 *
 * Lead developers: Max Sagebaum, Johannes Blühdorn (SciComp, University of Kaiserslautern-Landau)
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
 * For other licensing options please contact us.
 *
 * Authors:
 *  - SciComp, University of Kaiserslautern-Landau:
 *    - Max Sagebaum
 *    - Johannes Blühdorn
 *    - Former members:
 *      - Tim Albring
 */
#pragma once

#include <vector>

#include "../../config.h"
#include "../../misc/macros.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Ease of access structure for user-provided data on the tape for external functions. See
   * ExternalFunctionTapeInterface.
   *
   * Stores copies of the data provided to the add methods, either a single value or an entire array.
   *
   * The data can be retrieved in two different manners. Subsequent calls to get* methods provide the data elements in
   * the order in which they were added (In order access.). The get*ByIndex methods can be used to query the pos-th
   * added item explicitly (Out of order access.).
   *
   * After the last data element is accessed by a get* method, the next get* call will return the first one.
   *
   * The direct access function will not change the internal positioning.
   */
  struct ExternalFunctionUserData {
    private:

      struct DataItemBase {
        public:
          void* data;

          virtual ~DataItemBase() {}

          virtual DataItemBase* clone() = 0;
      };

      template<typename T_Type>
      struct DataItem : public DataItemBase {
        public:
          using Type = CODI_DD(T_Type, CODI_ANY);

          explicit DataItem(Type const& value) {
            data = (void*)new Type(value);
          }

          ~DataItem() {
            Type* pointer = (Type*)data;
            delete pointer;
          }

          DataItemBase* clone() {
            return new DataItem<Type>(*((Type*)data));
          }
      };

      template<typename Type>
      struct DataArray : public DataItemBase {
        private:
          int size;

        public:
          DataArray(Type const* value, int size) {
            data = (void*)new Type[size];
            this->size = size;
            std::copy(value, &value[size], (Type*)data);
          }

          ~DataArray() {
            Type* pointer = (Type*)data;
            delete[] pointer;
          }

          DataItemBase* clone() {
            return new DataArray<Type>(*((Type*)data), size);
          }
      };

      std::vector<DataItemBase*> store;

      size_t storePos;

    public:

      /// Constructor
      ExternalFunctionUserData() : storePos(0) {}

      /// Constructor
      ExternalFunctionUserData(ExternalFunctionUserData const& other) {
        copyAll(other, *this);
      }

      /// Destructor
      ~ExternalFunctionUserData() {
        clear();
      }

      /// Copy operator
      ExternalFunctionUserData& operator=(ExternalFunctionUserData const& other) {
        this->clear();
        copyAll(other, *this);

        return *this;
      }

      /// Delete all data entries.
      void clear() {
        for (size_t i = 0; i < store.size(); ++i) {
          delete store[i];
        }
        store.clear();
      }

      /// Add a value to the store. The value is copied.
      ///
      /// @return Index of the value for direct access.
      template<typename Type>
      size_t addData(Type const& value) {
        store.push_back(new DataItem<Type>(value));
        return store.size() - 1;
      }

      /// Add an array to the store. The array is copied.
      ///
      /// @return Index of the array for direct access.
      template<typename Type>
      size_t addData(Type const* value, int const size) {
        store.push_back(new DataArray<Type>(value, size));
        return store.size() - 1;
      }

      /*******************************************************************************/
      /// @name In order accessors
      /// @{

      /// Get a copy of the next data item.
      template<typename Type>
      void getData(Type& value) {
        getData<Type>(&value, 1);
      }

      /// Get a constant reference to the next data item.
      template<typename Type>
      Type const& getData() {
        Type* data = nextStore<Type>();
        return *data;
      }

      /// Get a reference to the next data item.
      template<typename Type>
      Type& getDataRef() {
        Type* data = nextStore<Type>();
        return *data;
      }

      /// Get the next data item and copy it as an array. The target array must have the correct size.
      template<typename Type>
      void getData(Type* value, int const size) {
        Type* convPointer = nextStore<Type>();

        std::copy(convPointer, &convPointer[size], value);
      }

      /// Get the address of the next data item. Intended for reference access to array items.
      template<typename Type>
      Type const* getDataArray() {
        Type* data = nextStore<Type>();
        return *data;
      }

      /// Manually reset the position.
      void resetPos() {
        storePos = 0;
      }

      /// @}

      /*******************************************************************************/
      /// @name Out of order accessors
      /// @{

      /// Get a copy of a data item based on the index.
      ///
      /// The internal position is not modified.
      template<typename Type>
      void getDataByIndex(Type& value, size_t pos) {
        getDataByIndex<Type>(&value, 1, pos);
      }

      /// Get a constant reference to a data item based on the index.
      ///
      /// The internal position is not modified.
      template<typename Type>
      Type const& getDataByIndex(size_t pos) {
        Type* data = getStore<Type>(pos);
        return *data;
      }

      /// Get a reference to a item based on the index.
      ///
      /// The internal position is not modified.
      template<typename Type>
      Type& getDataRefByIndex(size_t pos) {
        Type* data = getStore<Type>(pos);
        return *data;
      }

      /// Get a data item based on the index and copy it as an array. The target array must have the correct size.
      ///
      /// The internal position is not modified.
      template<typename Type>
      void getDataByIndex(Type* value, int const size, size_t pos) {
        Type* convPointer = getStore<Type>(pos);

        std::copy(convPointer, &convPointer[size], value);
      }

      /// Get the address of a data item based on the index. Intended for reference access to array items.
      template<typename Type>
      Type const* getDataArrayByIndex() {
        Type* data = nextStore<Type>();
        return *data;
      }

      /// @}

    private:

      static void copyAll(ExternalFunctionUserData const& from, ExternalFunctionUserData& to) {
        for (size_t i = 0; i < from.store.size(); ++i) {
          to.store.push_back(from.store[i]->clone());
        }
        to.storePos = from.storePos;
      }

      template<typename Type>
      Type* getStore(size_t pos) {
        return (Type*)store[pos]->data;
      }

      template<typename Type>
      Type* nextStore() {
        Type* pointer = getStore<Type>(storePos);
        storePos += 1;
        if (storePos >= store.size()) {
          storePos = 0;
        }

        return pointer;
      }
  };
}
