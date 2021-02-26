#pragma once

#include <vector>

#include "../../aux/macros.hpp"
#include "../../config.h"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Ease of access structure for used data on the tape for external functions. See ExternalFunctionTapeInterface.
   *
   * The structure clones all data that is stored via the put methods.
   *
   * Each get will retrieve the data in the same order as it was put into the structure. After the last data element is
   * accessed the first one will be returned on the next get.
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

      template<typename _Type>
      struct DataItem : public DataItemBase {
        public:
          using Type = CODI_DD(_Type, CODI_ANY);

          explicit DataItem(Type const& value) {
            data = (void*) new Type(value);
          }

          ~DataItem() {
            Type* pointer = (Type*) data;
            delete pointer;
          }

          DataItemBase* clone() { return new DataItem<Type>(*((Type*) data)); }
      };

      template<typename Type>
      struct DataArray : public DataItemBase {
        private:
          int size;
        public:
          DataArray(Type const* value, int size) {
            data = (void*) new Type[size];
            this->size = size;
            std::copy(value, &value[size], (Type*) data);
          }

          ~DataArray() {
            Type* pointer = (Type*) data;
            delete[] pointer;
          }

          DataItemBase* clone() { return new DataArray<Type>(*((Type*) data), size); }
      };

      std::vector<DataItemBase*> store;

      size_t storePos;

    public:

      /// Constructor
      ExternalFunctionUserData() :
        storePos(0) { }

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

      /// Delete all data entries
      void clear() {
        for (size_t i = 0; i < store.size(); ++i) {
          delete store[i];
        }
        store.clear();
      }

      /// Add a value to the store.
      ///
      /// @return index of the value for direct access
      template<typename Type>
      size_t addData(Type const& value) {
        store.push_back(new DataItem<Type>(value));
        return store.size() - 1;
      }

      /// Add an array to the store. The array is copied.
      ///
      /// @return index of the array for direct access
      template<typename Type>
      size_t addData(Type const* value, int const size) {
        store.push_back(new DataArray<Type>(value, size));
        return store.size() - 1;
      }

      /// Get the next data item and copy it into the value.
      template<typename Type>
      void getData(Type& value) {
        getData<Type>(&value, 1);
      }

      /// Get the next data item and return a constant reference
      template<typename Type>
      Type const& getData() {
        Type* data = nextStore<Type>();
        return *data;
      }

      /// Get a writable reference to the next data item.
      template<typename Type>
      Type& getDataRef() {
        Type* data = nextStore<Type>();
        return *data;
      }

      /// Get the next data item and copy the array. The target array must have the correct size.
      template<typename Type>
      void getData(Type* value, int const size) {
        Type* convPointer = nextStore<Type>();

        std::copy(convPointer, &convPointer[size], value);
      }

      /// The the next data item and return a constant reference
      template<typename Type>
      Type const* getDataArray() {
        Type* data = nextStore<Type>();
        return *data;
      }

      /// Get a data item based on the index.
      ///
      /// The internal position is not modified.
      template<typename Type>
      void getDataByIndex(Type& value, size_t pos) {
        getDataArrayByIndex<Type>(&value, 1, pos);
      }


      /// Get a data array based on the index. The target array must have the correct size.
      ///
      /// The internal position is not modified.
      template<typename Type>
      void getDataArrayByIndex(Type* value, int const size, size_t pos) {
        Type* convPointer = getStore<Type>(pos);

        std::copy(convPointer, &convPointer[size], value);
      }

      /// Manually reset the position
      void resetPos() {
        storePos = 0;
      }

    private:

      static void copyAll(ExternalFunctionUserData const& from, ExternalFunctionUserData& to) {
        for (size_t i = 0; i < from.store.size(); ++i) {
          to.store.push_back(from.store[i]->clone());
        }
        to.storePos = from.storePos;
      }

      template<typename Type>
      Type* getStore(size_t pos) {
        return (Type*) store[pos]->data;
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
