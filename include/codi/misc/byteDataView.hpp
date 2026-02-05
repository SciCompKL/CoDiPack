/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2026 Chair for Scientific Computing (SciComp), RPTU University Kaiserslautern-Landau
 * Homepage: http://scicomp.rptu.de
 * Contact:  Prof. Nicolas R. Gauger (codi@scicomp.uni-kl.de)
 *
 * Lead developers: Max Sagebaum, Johannes Blühdorn (SciComp, RPTU University Kaiserslautern-Landau)
 *
 * This file is part of CoDiPack (http://scicomp.rptu.de/software/codi).
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
 *  - SciComp, RPTU University Kaiserslautern-Landau:
 *    - Max Sagebaum
 *    - Johannes Blühdorn
 *    - Former members:
 *      - Tim Albring
 */
#pragma once

#include <algorithm>
#include <vector>

#include "../config.h"
#include "macros.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /// Helper for reading and writing to a byte array.
  ///
  /// Returned pointers are always pointers into the byte array. Every read, write, and reserve function advances the
  /// internal position in the byte array by the size of the data. It is always assumed that enough space is left in
  /// the byte array.
  struct ByteDataView {
    private:

      char* pointer;  ///< Data pointer.
      size_t pos;     ///< Current data position.

      size_t start;  ///< Start of the available data.
      size_t end;    ///< Size of the available data.

    public:

      /// Empty initialization. Innit needs to be called on this object.
      CODI_INLINE ByteDataView() = default;

      /// Constructor.
      CODI_INLINE ByteDataView(char* pointer, size_t pos, size_t end)
          : pointer(pointer), pos(pos), start(pos), end(end) {
        codiAssert(pos <= end);
      }

      /// Get the end data position.
      CODI_INLINE size_t getEnd() {
        return end;
      }

      /// Get the current data position.
      CODI_INLINE size_t getPosition() {
        return pos;
      }

      /// Get the start data position.
      CODI_INLINE size_t getStart() {
        return start;
      }

      /// Initialize the object.
      CODI_INLINE void init(char* pointer, size_t pos, size_t end) {
        this->pointer = pointer;
        this->pos = pos;
        this->start = pos;
        this->end = end;

        codiAssert(pos <= end);
      }

      /// Read an array of length \c size of types \c T.
      template<typename T>
      CODI_INLINE T* read(size_t size) {
        T* convPointer = cast<T>();
        pos += sizeof(T) * size;

        codiAssert(pos <= end);

        return convPointer;
      }

      /// Read a single object of type \c T.
      template<typename T>
      CODI_INLINE T read() {
        return read<T>(1)[0];
      }

      /// @brief Reserve memory for an array of type \c T with length \c size. The returned pointer can be written with
      /// the actual data.
      template<typename T>
      CODI_INLINE T* reserve(size_t size) {
        T* convPointer = cast<T>();
        pos += sizeof(T) * size;

        codiAssert(pos <= end);

        return convPointer;
      }

      /// @brief Reset the data position to the start of the data.
      CODI_INLINE void reset() {
        pos = start;
      }

      /// Write a single entry of type \c T.
      template<typename T>
      CODI_INLINE T* write(T const& data) {
        return write(&data, 1);
      }

      /// Write an array of type \c T with length \c size.
      template<typename T>
      CODI_INLINE T* write(T const* data, size_t size) {
        T* convPointer = cast<T>();
        for (size_t i = 0; i < size; i += 1) {
          convPointer[i] = data[i];
        }
        pos += sizeof(T) * size;

        codiAssert(pos <= end);

        return convPointer;
      }

    private:

      template<typename T>
      CODI_INLINE T* cast() {
        return reinterpret_cast<T*>(&pointer[pos]);
      }
  };

}
