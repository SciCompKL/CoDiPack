#pragma once

#include <vector>

#include "../../aux/macros.h"
#include "../../config.h"

/** \copydoc codi::Namespace */
namespace codi {

  struct ExternalFunctionData {
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
          using Type = DECLARE_DEFAULT(_Type, ANY);

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
          DataArray(const Type* value, int size) {
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

      ExternalFunctionData() :
        storePos(0) { }

      ~ExternalFunctionData() {
        clear();
      }

      void clear() {
        for(size_t i = 0; i < store.size(); ++i) {
          delete store[i];
        }
        store.clear();
      }

      ExternalFunctionData(ExternalFunctionData const& other) {
        copyAll(other, *this);
      }

      ExternalFunctionData& operator=(ExternalFunctionData const& other) {
        this->clear();
        copyAll(other, *this);

        return *this;
      }

      template<typename Type>
      size_t addData(Type const& value) {
        store.push_back(new DataItem<Type>(value));
        return store.size() - 1;
      }

      template<typename Type>
      size_t addData(Type const* value, int const size) {
        store.push_back(new DataArray<Type>(value, size));
        return store.size() - 1;
      }

      template<typename Type>
      void getData(Type& value) {
        getData<Type>(&value, 1);
      }

      template<typename Type>
      Type const& getData() {
        Type* data = nextStore<Type>();
        return *data;
      }

      template<typename Type>
      Type& getDataRef() {
        Type* data = nextStore<Type>();
        return *data;
      }

      template<typename Type>
      void getData(Type* value, int const size) {
        Type* convPointer = nextStore<Type>();

        std::copy(convPointer, &convPointer[size], value);
      }

      template<typename Type>
      Type const* getDataArray() {
        Type* data = nextStore<Type>();
        return *data;
      }

      template<typename Type>
      void getDataByIndex(Type& value, size_t pos) {
        getDataArrayByIndex<Type>(&value, 1, pos);
      }

      template<typename Type>
      void getDataArrayByIndex(Type* value, int const size, size_t pos) {
        Type* convPointer = getStore<Type>(pos);

        std::copy(convPointer, &convPointer[size], value);
      }

      void resetPos() {
        storePos = 0;
      }

    private:

      static void copyAll(ExternalFunctionData const& from, ExternalFunctionData& to) {
        for(size_t i = 0; i < from.store.size(); ++i) {
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
        if(storePos >= store.size()) {
          storePos = 0;
        }

        return pointer;
      }
  };
}
