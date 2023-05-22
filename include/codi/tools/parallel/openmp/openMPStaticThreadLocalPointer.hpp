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

#include "../staticThreadLocalPointerInterface.hpp"
#include "macros.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Static thread-local pointers for OpenMP.
   *
   * @tparam T_Type   See StaticThreadLocalPointerInterface.
   * @tparam T_Owner  See StaticThreadLocalPointerInterface.
   */
  template<typename T_Type, typename T_Owner>
  struct OpenMPStaticThreadLocalPointer
      : public StaticThreadLocalPointerInterface<T_Type, T_Owner, OpenMPStaticThreadLocalPointer<T_Type, T_Owner>> {
    public:
      using Type = T_Type;                       ///< See OpenMPStaticThreadLocalPointer.
      using Owner = CODI_DD(T_Owner, CODI_ANY);  ///< See OpenMPStaticThreadLocalPointer.

    private:
      // Returning the static thread local pointer from a function works around a tls bug of gcc,
      // see https://gcc.gnu.org/bugzilla/show_bug.cgi?id=66944 for details.
      static CODI_INLINE Type*& getPtr() {
        static Type* value = new Type();
        CODI_OMP_THREADPRIVATE(value)

        return value;
      }

    public:

      /// \copydoc StaticThreadLocalPointerInterface::set
      static CODI_INLINE void set(Type* other) {
        getPtr() = other;
      }

      /// \copydoc StaticThreadLocalPointerInterface::get
      static CODI_INLINE Type* get() {
        return getPtr();
      }
  };
}
