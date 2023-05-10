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

#include <cstddef>

#include "../../config.h"
#include "../../misc/fileIo.hpp"
#include "../../misc/macros.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief A chunk stores a contiguous block of data in CoDiPack.
   *
   * See DataInterface for a more general description of the data layout in CoDiPack.
   *
   * The chunk classes implement a structure of arrays approach for the data management. Each item can have multiple
   * entries where each entry is stored in its own array.
   *
   * E.g., if each item consists two entries (double, int), then we have two arrays:
   *
   * \code{.txt}
   *                     item 0 | item 1 | item 2 | etc.
   *   array1 (double) :  0.1   |   3.14 |  2.17  | ...
   *   array2 (int)    :   1    |   10   |   2    | ...
   * \endcode
   *
   * ChunkBase serves both as a base class and as an interface.
   * The base class defines functions for getting and setting the number of used items. The interface defines the
   * functions for the data access.
   *
   * - Entry management:
   *   - pushData(): Add one data item. One argument per entry.
   *   - dataPointer(): Get pointers to the data. One argument per entry.
   *   - erase(): Erase a range of data items.
   *
   * - Data IO:
   *   - allocateData() / deleteData(): Allocate / delete the data arrays.
   *   - readData() / writeData(): Read / write the data in the arrays to the IO object.
   *
   */
  struct ChunkBase {
    public:

      /*******************************************************************************/
      /// @name Interface: Types & constants
      /// @{

      static size_t constexpr EntrySize = CODI_UNDEFINED_VALUE;  ///< Total size of all data in one entry.

      /// @}
      /*******************************************************************************/
      /// @name Interface: Entry management
      /// @{

      template<typename... Data>
      CODI_INLINE void pushData(Data&&... dataEntries);  ///< Add one data item. For each entry one argument has to be
                                                         ///< provided.

      template<typename... Pointers>
      CODI_INLINE void dataPointer(size_t const& index,
                                   Pointers*&... pointers);  ///< Extract pointer to requested position. For each entry
                                                             ///< one argument has to be provided.

      virtual void erase(size_t const& start, size_t const& end) = 0;  ///< Erase data items start, ..., end - 1.

      /// @}
      /*******************************************************************************/
      /// @name Interface: Data IO
      /// @{

      virtual void allocateData() = 0;                   ///< Allocated the data if it was deallocated before.
      virtual void deleteData() = 0;                     ///< Delete the allocated data.
      virtual void readData(FileIo& handle) = 0;         ///< Read data from the FileIo handle.
      virtual void writeData(FileIo& handle) const = 0;  ///< Write data to the FileIo handle.

      /// @}
      /*******************************************************************************/
      /// @name Interface: Misc
      /// @{

      void swap(CODI_IMPLEMENTATION& other);  ///< Swap data with other chunk of the same type.

      /// @}
    protected:
      size_t size;      ///< Maximum size of arrays.
      size_t usedSize;  ///< Currently used size.

    public:

      /// Constructor
      explicit ChunkBase(size_t const& size) : size(size), usedSize(0) {}

      /// Destructor
      virtual ~ChunkBase() {}

      /*******************************************************************************/
      /// @name Common functionality
      /// @{

      /// Get the allocated size.
      CODI_INLINE size_t getSize() const {
        return size;
      }

      /// Number of unused data items.
      CODI_INLINE size_t getUnusedSize() const {
        return size - usedSize;
      }

      /// Number of used data items.
      CODI_INLINE size_t getUsedSize() const {
        return usedSize;
      }

      /// Sets the number of used items to zero.
      CODI_INLINE void reset() {
        usedSize = 0;
      }

      /// Resize the allocated data. Stored data is lost. Used size is set to zero.
      CODI_INLINE void resize(size_t newSize) {
        deleteData();
        size = newSize;
        usedSize = 0;
        allocateData();
      }

      /// Set the used size.
      CODI_INLINE void setUsedSize(size_t const& usage) {
        usedSize = usage;
      }

      /// @}

    protected:

      /// Swap the entries of this base class.
      void swap(ChunkBase& other) {
        std::swap(size, other.size);
        std::swap(usedSize, other.usedSize);
      }
  };

  /**
   * @ brief Chunk with one entry per item.
   *
   * @tparam Chunk1  Any type.
   */
  template<typename Data1>
  struct Chunk1 final : public ChunkBase {
    public:

      using Base = ChunkBase;  ///< Abbreviation for the base class type.

    private:

      Data1* data1;

    public:

      /// Constructor
      Chunk1(size_t const& size) : ChunkBase(size), data1(nullptr) {
        allocateData();
      }

      /// Destructor
      ~Chunk1() {
        deleteData();
      }

      /*******************************************************************************/
      /// @name ChunkBase interface implementation
      /// @{

      static size_t constexpr EntrySize = sizeof(Data1);  ///< \copydoc ChunkBase::EntrySize

      /// \copydoc ChunkBase::allocateData()
      void allocateData() {
        if (nullptr == data1) {
          data1 = new Data1[size];
        }
      }

      /// \copydoc ChunkBase::dataPointer
      CODI_INLINE void dataPointer(size_t const& index, Data1*& pointer1) {
        codiAssert(index <= ChunkBase::size);
        pointer1 = &data1[index];
      }

      /// \copydoc ChunkBase::deleteData
      void deleteData() {
        if (nullptr != data1) {
          delete[] data1;
          data1 = nullptr;
        }
      }

      /// \copydoc ChunkBase::erase
      void erase(size_t const& start, size_t const& end) {
        codiAssert(start <= end);
        codiAssert(start < usedSize);
        codiAssert(end <= usedSize);

        if (start != end) {
          for (size_t i = 0; i < usedSize - end; ++i) {
            data1[start + i] = data1[end + i];
          }
          usedSize -= end - start;
        }
      }

      /// \copydoc ChunkBase::pushData
      CODI_INLINE void pushData(Data1 const& value1) {
        codiAssert(getUnusedSize() != 0);
        data1[usedSize] = value1;
        usedSize += 1;
      }

      /// \copydoc ChunkBase::readData
      void readData(FileIo& handle) {
        allocateData();

        handle.readData(data1, size);
      }

      /// \copydoc ChunkBase::swap
      void swap(Chunk1<Data1>& other) {
        Base::swap(other);

        std::swap(data1, other.data1);
      }

      /// \copydoc ChunkBase::writeData
      void writeData(FileIo& handle) const {
        handle.writeData(data1, size);
      }

      /// @}
  };

  /**
   * @ brief Chunk with two entries per item.
   *
   * @tparam Chunk1  Any type.
   * @tparam Chunk2  Any type.
   */
  template<typename Data1, typename Data2>
  struct Chunk2 final : public ChunkBase {
    public:

      using Base = ChunkBase;  ///< Abbreviation for the base class type.

    private:

      Data1* data1;
      Data2* data2;

    public:

      /// Constructor
      Chunk2(size_t const& size) : ChunkBase(size), data1(nullptr), data2(nullptr) {
        allocateData();
      }

      /// Destructor
      ~Chunk2() {
        deleteData();
      }

      /*******************************************************************************/
      /// @name ChunkBase interface implementation
      /// @{

      static size_t constexpr EntrySize = sizeof(Data1) + sizeof(Data2);  ///< \copydoc ChunkBase::EntrySize

      /// \copydoc ChunkBase::allocateData()
      void allocateData() {
        if (nullptr == data1) {
          data1 = new Data1[size];
        }

        if (nullptr == data2) {
          data2 = new Data2[size];
        }
      }

      /// \copydoc ChunkBase::dataPointer
      CODI_INLINE void dataPointer(size_t const& index, Data1*& pointer1, Data2*& pointer2) {
        codiAssert(index <= ChunkBase::size);
        pointer1 = &data1[index];
        pointer2 = &data2[index];
      }

      /// \copydoc ChunkBase::deleteData
      void deleteData() {
        if (nullptr != data1) {
          delete[] data1;
          data1 = nullptr;
        }

        if (nullptr != data2) {
          delete[] data2;
          data2 = nullptr;
        }
      }

      /// \copydoc ChunkBase::erase
      void erase(size_t const& start, size_t const& end) {
        codiAssert(start <= end);
        codiAssert(start < usedSize);
        codiAssert(end <= usedSize);

        if (start != end) {
          for (size_t i = 0; i < usedSize - end; ++i) {
            data1[start + i] = data1[end + i];
          }
          for (size_t i = 0; i < usedSize - end; ++i) {
            data2[start + i] = data2[end + i];
          }
          usedSize -= end - start;
        }
      }

      /// \copydoc ChunkBase::pushData
      CODI_INLINE void pushData(Data1 const& value1, Data2 const& value2) {
        codiAssert(getUnusedSize() != 0);
        data1[usedSize] = value1;
        data2[usedSize] = value2;
        usedSize += 1;
      }

      /// \copydoc ChunkBase::readData
      void readData(FileIo& handle) {
        allocateData();

        handle.readData(data1, size);
        handle.readData(data2, size);
      }

      /// \copydoc ChunkBase::swap
      void swap(Chunk2<Data1, Data2>& other) {
        Base::swap(other);

        std::swap(data1, other.data1);
        std::swap(data2, other.data2);
      }

      /// \copydoc ChunkBase::writeData
      void writeData(FileIo& handle) const {
        handle.writeData(data1, size);
        handle.writeData(data2, size);
      }

      /// @}
  };

  /**
   * @ brief Chunk with three entries per item.
   *
   * @tparam Chunk1  Any type.
   * @tparam Chunk2  Any type.
   * @tparam Chunk3  Any type.
   */
  template<typename Data1, typename Data2, typename Data3>
  struct Chunk3 final : public ChunkBase {
    public:

      using Base = ChunkBase;  ///< Abbreviation for the base class type.

    private:

      Data1* data1;
      Data2* data2;
      Data3* data3;

    public:

      /// Constructor
      Chunk3(size_t const& size) : ChunkBase(size), data1(nullptr), data2(nullptr), data3(nullptr) {
        allocateData();
      }

      /// Destructor
      ~Chunk3() {
        deleteData();
      }

      /*******************************************************************************/
      /// @name ChunkBase interface implementation
      /// @{

      static size_t constexpr EntrySize =
          sizeof(Data1) + sizeof(Data2) + sizeof(Data3);  ///< \copydoc ChunkBase::EntrySize

      /// \copydoc ChunkBase::allocateData()
      void allocateData() {
        if (nullptr == data1) {
          data1 = new Data1[size];
        }

        if (nullptr == data2) {
          data2 = new Data2[size];
        }

        if (nullptr == data3) {
          data3 = new Data3[size];
        }
      }

      /// \copydoc ChunkBase::dataPointer
      CODI_INLINE void dataPointer(size_t const& index, Data1*& pointer1, Data2*& pointer2, Data3*& pointer3) {
        codiAssert(index <= ChunkBase::size);
        pointer1 = &data1[index];
        pointer2 = &data2[index];
        pointer3 = &data3[index];
      }

      /// \copydoc ChunkBase::deleteData
      void deleteData() {
        if (nullptr != data1) {
          delete[] data1;
          data1 = nullptr;
        }

        if (nullptr != data2) {
          delete[] data2;
          data2 = nullptr;
        }

        if (nullptr != data3) {
          delete[] data3;
          data3 = nullptr;
        }
      }

      /// \copydoc ChunkBase::erase
      void erase(size_t const& start, size_t const& end) {
        codiAssert(start <= end);
        codiAssert(start < usedSize);
        codiAssert(end <= usedSize);

        if (start != end) {
          for (size_t i = 0; i < usedSize - end; ++i) {
            data1[start + i] = data1[end + i];
          }
          for (size_t i = 0; i < usedSize - end; ++i) {
            data2[start + i] = data2[end + i];
          }
          for (size_t i = 0; i < usedSize - end; ++i) {
            data3[start + i] = data3[end + i];
          }
          usedSize -= end - start;
        }
      }

      /// \copydoc ChunkBase::pushData
      CODI_INLINE void pushData(Data1 const& value1, Data2 const& value2, Data3 const& value3) {
        codiAssert(getUnusedSize() != 0);
        data1[usedSize] = value1;
        data2[usedSize] = value2;
        data3[usedSize] = value3;
        usedSize += 1;
      }

      /// \copydoc ChunkBase::readData
      void readData(FileIo& handle) {
        allocateData();

        handle.readData(data1, size);
        handle.readData(data2, size);
        handle.readData(data3, size);
      }

      /// \copydoc ChunkBase::swap
      void swap(Chunk3<Data1, Data2, Data3>& other) {
        Base::swap(other);

        std::swap(data1, other.data1);
        std::swap(data2, other.data2);
        std::swap(data3, other.data3);
      }

      /// \copydoc ChunkBase::writeData
      void writeData(FileIo& handle) const {
        handle.writeData(data1, size);
        handle.writeData(data2, size);
        handle.writeData(data3, size);
      }

      /// @}
  };

  /**
   * @ brief Chunk with four entries per item.
   *
   * @tparam Chunk1  Any type.
   * @tparam Chunk2  Any type.
   * @tparam Chunk3  Any type.
   * @tparam Chunk4  Any type.
   */
  template<typename Data1, typename Data2, typename Data3, typename Data4>
  struct Chunk4 final : public ChunkBase {
    public:

      using Base = ChunkBase;  ///< Abbreviation for the base class type.

    private:

      Data1* data1;
      Data2* data2;
      Data3* data3;
      Data4* data4;

    public:

      /// Constructor
      Chunk4(size_t const& size) : ChunkBase(size), data1(nullptr), data2(nullptr), data3(nullptr), data4(nullptr) {
        allocateData();
      }

      /// Destructor
      ~Chunk4() {
        deleteData();
      }

      /*******************************************************************************/
      /// @name ChunkBase interface implementation
      /// @{

      static size_t constexpr EntrySize =
          sizeof(Data1) + sizeof(Data2) + sizeof(Data3) + sizeof(Data4);  ///< \copydoc ChunkBase::EntrySize

      /// \copydoc ChunkBase::allocateData()
      void allocateData() {
        if (nullptr == data1) {
          data1 = new Data1[size];
        }

        if (nullptr == data2) {
          data2 = new Data2[size];
        }

        if (nullptr == data3) {
          data3 = new Data3[size];
        }

        if (nullptr == data4) {
          data4 = new Data4[size];
        }
      }

      /// \copydoc ChunkBase::dataPointer
      CODI_INLINE void dataPointer(size_t const& index, Data1*& pointer1, Data2*& pointer2, Data3*& pointer3,
                                   Data4*& pointer4) {
        codiAssert(index <= ChunkBase::size);
        pointer1 = &data1[index];
        pointer2 = &data2[index];
        pointer3 = &data3[index];
        pointer4 = &data4[index];
      }

      /// \copydoc ChunkBase::deleteData
      void deleteData() {
        if (nullptr != data1) {
          delete[] data1;
          data1 = nullptr;
        }

        if (nullptr != data2) {
          delete[] data2;
          data2 = nullptr;
        }

        if (nullptr != data3) {
          delete[] data3;
          data3 = nullptr;
        }

        if (nullptr != data4) {
          delete[] data4;
          data4 = nullptr;
        }
      }

      /// \copydoc ChunkBase::erase
      void erase(size_t const& start, size_t const& end) {
        codiAssert(start <= end);
        codiAssert(start < usedSize);
        codiAssert(end <= usedSize);

        if (start != end) {
          for (size_t i = 0; i < usedSize - end; ++i) {
            data1[start + i] = data1[end + i];
          }
          for (size_t i = 0; i < usedSize - end; ++i) {
            data2[start + i] = data2[end + i];
          }
          for (size_t i = 0; i < usedSize - end; ++i) {
            data3[start + i] = data3[end + i];
          }
          for (size_t i = 0; i < usedSize - end; ++i) {
            data4[start + i] = data4[end + i];
          }
          usedSize -= end - start;
        }
      }

      /// \copydoc ChunkBase::pushData
      CODI_INLINE void pushData(Data1 const& value1, Data2 const& value2, Data3 const& value3, Data4 const& value4) {
        codiAssert(getUnusedSize() != 0);
        data1[usedSize] = value1;
        data2[usedSize] = value2;
        data3[usedSize] = value3;
        data4[usedSize] = value4;
        usedSize += 1;
      }

      /// \copydoc ChunkBase::readData
      void readData(FileIo& handle) {
        allocateData();

        handle.readData(data1, size);
        handle.readData(data2, size);
        handle.readData(data3, size);
        handle.readData(data4, size);
      }

      /// \copydoc ChunkBase::swap
      void swap(Chunk4<Data1, Data2, Data3, Data4>& other) {
        Base::swap(other);

        std::swap(data1, other.data1);
        std::swap(data2, other.data2);
        std::swap(data3, other.data3);
        std::swap(data4, other.data4);
      }

      /// \copydoc ChunkBase::writeData
      void writeData(FileIo& handle) const {
        handle.writeData(data1, size);
        handle.writeData(data2, size);
        handle.writeData(data3, size);
        handle.writeData(data4, size);
      }

      /// @}
  };

}
