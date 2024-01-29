/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2024 Chair for Scientific Computing (SciComp), University of Kaiserslautern-Landau
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

#include "../../../config.h"
#include "../../../misc/byteDataView.hpp"
#include "../../../misc/macros.hpp"
#include "../../../misc/temporaryMemory.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   *  @brief Traits for storing passive arguments in byte streams.
   *
   *  See ActiveArgumentStoreTraits for a detailed documentation of the process.
   *
   *  If \c storeRequired is false, nothing needs to be stored for this type and only a default initialization should
   *  be performed for the restore operation.
   *
   *  @tparam T_T  Type of the argument that is stored.
   *  @tparam T_S  Store type of the argument that is stored.
   */
  template<typename T_T, typename T_S = T_T, typename = void>
  struct PassiveArgumentStoreTraits {
      using T = CODI_DD(T_T, int);    ///< See PassiveArgumentStoreTraits.
      using S = CODI_DD(T_S, short);  ///< See PassiveArgumentStoreTraits.

      using Store = T;  ///< Type for the variable declaration for restoring the data.

      /// Count the required size for storing the data.
      CODI_INLINE static size_t countSize(T&& value, size_t size, bool storeRequired);

      /// Restore the data for this type.
      CODI_INLINE static void restore(ByteDataView* store, TemporaryMemory& allocator, size_t size, bool storeRequired,
                                      Store& value);

      /// Store the data for the type in the data store.
      CODI_INLINE static void store(ByteDataView* dataStore, TemporaryMemory& allocator, T const& value, size_t size,
                                    bool storeRequired);
  };

#ifndef DOXYGEN_DISABLE

  /// Specialization of PassiveArgumentStoreTraits for integral values.
  template<typename T_T, typename T_S>
  struct PassiveArgumentStoreTraits<T_T, T_S, typename std::enable_if<std::is_integral<T_T>::value>::type> {
      using T = CODI_DD(T_T, int);    ///< See PassiveArgumentStoreTraits.
      using S = CODI_DD(T_S, short);  ///< See PassiveArgumentStoreTraits.

      using Store = T;

      /// @copydoc PassiveArgumentStoreTraits::countSize()
      CODI_INLINE static size_t countSize(T const& value, size_t size, bool storeRequired) {
        CODI_UNUSED(value, size);

        size_t storeSize = 0;

        if (storeRequired) {
          storeSize += sizeof(S);
        }

        return storeSize;
      }

      /// @copydoc PassiveArgumentStoreTraits::restore()
      CODI_INLINE static void restore(ByteDataView* dataStore, TemporaryMemory& allocator, size_t size,
                                      bool storeRequired, Store& value) {
        CODI_UNUSED(allocator, size);
        if (storeRequired) {
          value = *dataStore->template read<S>(1);
        }
      }

      /// @copydoc PassiveArgumentStoreTraits::store()
      CODI_INLINE static void store(ByteDataView* dataStore, TemporaryMemory& allocator, T const& value, size_t size,
                                    bool storeRequired) {
        CODI_UNUSED(allocator, size);

        if (storeRequired) {
          dataStore->write<S>(value);
        }
      }
  };

#endif
}
