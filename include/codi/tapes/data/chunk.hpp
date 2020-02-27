#pragma once

#include <cstddef>

#include "../../aux/macros.h"
#include "../../config.h"

/** \copydoc codi::Namespace */
namespace codi {

  struct ChunkBase {
    public:

      /*******************************************************************************
       * Section: Definition of the interface
       *
       * Description: TODO
       *
       */

      static size_t constexpr EntrySize = UNDEFINED_VALUE;

      void swap(ChunkBase& other);

      template<typename ... Data>
      CODI_INLINE void pushData(Data&& ... dataEntries);

      template<typename ... Pointers>
      CODI_INLINE void dataPointer(size_t const& index, Pointers*& ... pointers);

      //TODO: Add IO interface

      /*******************************************************************************
       * Section: Implementation of common functionality
       *
       * Description: TODO
       *
       */

      size_t size;
      size_t usedSize;

      explicit ChunkBase(size_t const& size) :
        size(size),
        usedSize(0) {}

      virtual ~ChunkBase() {}


      CODI_INLINE size_t getSize() const {
        return size;
      }

      CODI_INLINE size_t getUnusedSize() const {
        return size - usedSize;
      }

      CODI_INLINE size_t getUsedSize() const {
        return usedSize;
      }

      CODI_INLINE void reset() {
        usedSize = 0;
      }

      CODI_INLINE void setUsedSize(size_t const& usage) {
        usedSize = usage;
      }

    protected:

      void swapBase(ChunkBase& other) {
        std::swap(size, other.size);
        std::swap(usedSize, other.usedSize);
      }
  };

  template<typename Data1>
  struct Chunk1 final : public ChunkBase {
    public:

      static size_t constexpr EntrySize = sizeof(Data1);

      Data1* data1;

      Chunk1(size_t const& size) : ChunkBase(size),
        data1(NULL) {

        allocateData();
      }

      ~Chunk1() {
        deleteData();
      }

      void allocateData() {
        if(NULL == data1) {
          data1 = new Data1[size];
        }
      }

      CODI_INLINE void dataPointer(size_t const& index, Data1* &pointer1) {
        codiAssert(index <= ChunkBase::size);
        pointer1 = &data1[index];
      }

      void deleteData() {
        if(NULL != data1) {
          delete [] data1;
          data1 = NULL;
        }
      }

      CODI_INLINE void pushData(Data1 const& value1) {
        codiAssert(getUnusedSize() != 0);
        data1[usedSize] = value1;
        usedSize += 1;
      }

      void swap(Chunk4<Data1>& other) {
        this->swapBase(other);

        std::swap(data1, other.data1);
      }
  };

  template<typename Data1, typename Data2>
  struct Chunk2 final : public ChunkBase {
    public:

      static size_t constexpr EntrySize = sizeof(Data1) + sizeof(Data2);

      Data1* data1;
      Data2* data2;

      Chunk2(size_t const& size) : ChunkBase(size),
        data1(NULL),
        data2(NULL) {

        allocateData();
      }

      ~Chunk2() {
        deleteData();
      }

      void allocateData() {
        if(NULL == data1) {
          data1 = new Data1[size];
        }

        if(NULL == data2) {
          data2 = new Data2[size];
        }
      }

      CODI_INLINE void dataPointer(size_t const& index, Data1* &pointer1, Data2* &pointer2) {
        codiAssert(index <= ChunkBase::size);
        pointer1 = &data1[index];
        pointer2 = &data2[index];
      }

      void deleteData() {
        if(NULL != data1) {
          delete [] data1;
          data1 = NULL;
        }

        if(NULL != data2) {
          delete [] data2;
          data2 = NULL;
        }
      }

      CODI_INLINE void pushData(Data1 const& value1, Data2 const& value2) {
        codiAssert(getUnusedSize() != 0);
        data1[usedSize] = value1;
        data2[usedSize] = value2;
        usedSize += 1;
      }

      void swap(Chunk4<Data1, Data2>& other) {
        this->swapBase(other);

        std::swap(data1, other.data1);
        std::swap(data2, other.data2);
      }
  };

  template<typename Data1, typename Data2, typename Data3>
  struct Chunk3 final : public ChunkBase {
    public:

      static size_t constexpr EntrySize = sizeof(Data1) + sizeof(Data2) + sizeof(Data3);

      Data1* data1;
      Data2* data2;
      Data3* data3;

      Chunk3(size_t const& size) : ChunkBase(size),
        data1(NULL),
        data2(NULL),
        data3(NULL) {

        allocateData();
      }

      ~Chunk3() {
        deleteData();
      }

      void allocateData() {
        if(NULL == data1) {
          data1 = new Data1[size];
        }

        if(NULL == data2) {
          data2 = new Data2[size];
        }

        if(NULL == data3) {
          data3 = new Data3[size];
        }
      }

      CODI_INLINE void dataPointer(size_t const& index, Data1* &pointer1, Data2* &pointer2, Data3* &pointer3) {
        codiAssert(index <= ChunkBase::size);
        pointer1 = &data1[index];
        pointer2 = &data2[index];
        pointer3 = &data3[index];
      }

      void deleteData() {
        if(NULL != data1) {
          delete [] data1;
          data1 = NULL;
        }

        if(NULL != data2) {
          delete [] data2;
          data2 = NULL;
        }

        if(NULL != data3) {
          delete [] data3;
          data3 = NULL;
        }
      }

      CODI_INLINE void pushData(Data1 const& value1, Data2 const& value2, Data3 const& value3) {
        codiAssert(getUnusedSize() != 0);
        data1[usedSize] = value1;
        data2[usedSize] = value2;
        data3[usedSize] = value3;
        usedSize += 1;
      }

      void swap(Chunk4<Data1, Data2, Data3>& other) {
        this->swapBase(other);

        std::swap(data1, other.data1);
        std::swap(data2, other.data2);
        std::swap(data3, other.data3);
      }
  };


  template<typename Data1, typename Data2, typename Data3, typename Data4>
  struct Chunk4 final : public ChunkBase {
    public:

      static size_t constexpr EntrySize = sizeof(Data1) + sizeof(Data2) + sizeof(Data3) + sizeof(Data4);

      Data1* data1;
      Data2* data2;
      Data3* data3;
      Data4* data4;

      Chunk4(size_t const& size) : ChunkBase(size),
        data1(NULL),
        data2(NULL),
        data3(NULL),
        data4(NULL) {

        allocateData();
      }

      ~Chunk4() {
        deleteData();
      }

      void allocateData() {
        if(NULL == data1) {
          data1 = new Data1[size];
        }

        if(NULL == data2) {
          data2 = new Data2[size];
        }

        if(NULL == data3) {
          data3 = new Data3[size];
        }

        if(NULL == data4) {
          data4 = new Data4[size];
        }
      }

      CODI_INLINE void dataPointer(size_t const& index, Data1* &pointer1, Data2* &pointer2, Data3* &pointer3, Data4* &pointer4) {
        codiAssert(index <= ChunkBase::size);
        pointer1 = &data1[index];
        pointer2 = &data2[index];
        pointer3 = &data3[index];
        pointer4 = &data4[index];
      }

      void deleteData() {
        if(NULL != data1) {
          delete [] data1;
          data1 = NULL;
        }

        if(NULL != data2) {
          delete [] data2;
          data2 = NULL;
        }

        if(NULL != data3) {
          delete [] data3;
          data3 = NULL;
        }

        if(NULL != data4) {
          delete [] data4;
          data4 = NULL;
        }
      }

      CODI_INLINE void pushData(Data1 const& value1, Data2 const& value2, Data3 const& value3, Data4 const& value4) {
        codiAssert(getUnusedSize() != 0);
        data1[usedSize] = value1;
        data2[usedSize] = value2;
        data3[usedSize] = value3;
        data4[usedSize] = value4;
        usedSize += 1;
      }

      void swap(Chunk4<Data1, Data2, Data3, Data4>& other) {
        this->swapBase(other);

        std::swap(data1, other.data1);
        std::swap(data2, other.data2);
        std::swap(data3, other.data3);
        std::swap(data4, other.data4);
      }
  };

}
