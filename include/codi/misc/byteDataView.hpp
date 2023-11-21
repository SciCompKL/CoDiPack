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

#include <algorithm>
#include <vector>

#include "../config.h"
#include "macros.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /// Helper for reading and writing to a byte array.
  ///
  /// Returned pointers are always pointers into the byte array. Every read, write and reserve function advances the
  /// position into the byte array in the correct direction and by the size of the data. It is always assumed that
  /// enough space is left in the byte array.
  struct ByteDataView {
      /// Direction for reading and writing.
      enum class Direction {
        Forward,
        Reverse
      };

    private:

      char* pointer;  ///< Data pointer.
      size_t pos;     ///< Current data position.

      Direction direction;  ///< Direction for reading.

    public:

      /// Empty initialization. Innit needs to be called on this object.
      CODI_INLINE ByteDataView() = default;

      /// Constructor.
      CODI_INLINE ByteDataView(char* pointer, size_t pos, Direction d) : pointer(pointer), pos(pos), direction(d) {}

      /// Get the current data position.
      CODI_INLINE size_t getPosition() {
        return pos;
      }

      /// Get the reading direction.
      CODI_INLINE Direction getDirection() {
        return direction;
      }

      /// Initialize the object.
      CODI_INLINE void init(char* pointer, size_t pos, Direction d) {
        this->pointer = pointer;
        this->pos = pos;
        this->direction = d;
      }

      /// Read an array of length \c size of types \c T.
      template<typename T>
      CODI_INLINE T* read(size_t size) {
        if (Direction::Reverse == direction) {
          pos -= sizeof(T) * size;
        }
        T* convPointer = cast<T>();
        if (Direction::Forward == direction) {
          pos += sizeof(T) * size;
        }

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

        return convPointer;
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

        return convPointer;
      }

    private:

      template<typename T>
      CODI_INLINE T* cast() {
        return reinterpret_cast<T*>(&pointer[pos]);
      }
  };

}
