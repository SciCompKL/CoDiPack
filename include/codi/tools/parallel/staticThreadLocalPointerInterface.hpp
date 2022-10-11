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

#include "../../misc/macros.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Abstracts a static threadlocal pointer.
   *
   * @tparam T_Type   The datatype pointed to.
   * @tparam T_Owner  Type that owns the static pointer. Needed to distinguish between multiple such pointers.
   * @tparam T_Impl   The implementing type.
   */
  template<typename T_Type, typename T_Owner, typename T_Impl>
  struct StaticThreadLocalPointerInterface {
    public:
      using Type = CODI_DD(T_Type, CODI_ANY);             ///< See StaticThreadLocalPointerInterface.
      using Owner = CODI_DD(T_Owner, CODI_ANY);           ///< See StaticThreadLocalPointerInterface.
      using Impl = CODI_DD(T_Impl, CODI_IMPLEMENTATION);  ///< See StaticThreadLocalPointerInterface.

      CODI_INLINE StaticThreadLocalPointerInterface() {}  ///< Constructor.
      ~StaticThreadLocalPointerInterface() {}             ///< Destructor.

      static CODI_INLINE void set(Type* other);  ///< Set the pointer.
      static CODI_INLINE Type* get();            ///< Get the pointer.´
  };
}
