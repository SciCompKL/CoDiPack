/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2022 Chair for Scientific Computing (SciComp), TU Kaiserslautern
 * Homepage: http://www.scicomp.uni-kl.de
 * Contact:  Prof. Nicolas R. Gauger (codi@scicomp.uni-kl.de)
 *
 * Lead developers: Max Sagebaum, Johannes Blühdorn (SciComp, TU Kaiserslautern)
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
 *  - SciComp, TU Kaiserslautern:
 *    - Max Sagebaum
 *    - Johannes Blühdorn
 *    - Former members:
 *      - Tim Albring
 */
#pragma once

#include "../staticThreadLocalPointerInterface.hpp"

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
      using Type = CODI_DD(T_Type, CODI_ANY);    ///< See OpenMPStaticThreadLocalPointer.
      using Owner = CODI_DD(T_Owner, CODI_ANY);  ///< See OpenMPStaticThreadLocalPointer.

    private:
      static Type* value;
      #pragma omp threadprivate(value)

    public:

      /// \copydoc StaticThreadLocalPointerInterface::set
      static CODI_INLINE void set(Type* other) {
        value = other;
      }

      /// \copydoc StaticThreadLocalPointerInterface::get
      static CODI_INLINE Type* get() {
        return value;
      }
  };

  template<typename Type, typename Owner>
  Type* OpenMPStaticThreadLocalPointer<Type, Owner>::value = new Type();
}
