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

  /**
   *  @brief Allocator for temporary used memory.
   *
   *  Can be used in places where memory is often allocated and deallocated. This reduces reduce the overhead of the
   *  system calls.
   *
   *  The initial memory is 4 MiB and can be extended with a call to ensureSize. All memory is initialized with zeros.
   */
  struct TemporaryMemory {
      static size_t constexpr InitialDataSize = 4 * 1024 * 1024;  ///< 4 MiB of memory.

    private:

      std::vector<char> data;  ///< Allocated data.
      size_t dataPos;          ///< Current data position.

    public:

      /// Constructor.
      CODI_INLINE TemporaryMemory() : data(InitialDataSize), dataPos() {}

      /// Constructor.
      CODI_INLINE TemporaryMemory(size_t initialSize) : data(initialSize), dataPos() {}

      /// Returns true if no data is currently allocated.
      CODI_INLINE bool isEmpty() {
        return 0 == dataPos;
      }

      /// @brief Allocate an array of type \c T with length \c size. Data is zero initialized. No constructors of \c T
      /// are called.
      template<typename T>
      CODI_INLINE T* alloc(size_t size) {
        codiAssert(dataPos + size * sizeof(T) <= data.size());

        T* castPointer = reinterpret_cast<T*>(&data.data()[dataPos]);
        dataPos += size * sizeof(T);

        return castPointer;
      }

      /// Allocate a single entity of \c T and call the constructor with \c args.
      template<typename T, typename... Args>
      CODI_INLINE T* allocAndInit(Args&&... args) {
        T* value = alloc<T>(1);
        new (value) T(std::forward<Args>(args)...);

        return value;
      }

      /// Ensures that enough space is available. Can only be called when no data has been allocated because
      /// reallocations invalidate pointers.
      CODI_INLINE void ensureSize(size_t newSize) {
        if (dataPos != 0) {
          CODI_EXCEPTION("Temporary memory can only be extended when no data is allocated.");
        }

        if (data.size() < newSize) {
          data.resize(newSize);
        }
      }

      /// @brief Free all allocated memory. No destructors are called. Stored pointers and resources need to be
      /// deallocated manually beforehand.
      CODI_INLINE void free() {
        // Clear used data.
        std::fill(data.begin(), data.begin() + dataPos, 0);
        dataPos = 0;
      }
  };
}
