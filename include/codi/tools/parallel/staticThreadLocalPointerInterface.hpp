/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2021 Chair for Scientific Computing (SciComp), TU Kaiserslautern
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

#include "../../../misc/macros.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  template<typename T_Type, typename T_Owner, typename T_Impl>
  struct StaticThreadLocalPointerInterface {
    public:
      using Type = CODI_DD(T_Type, CODI_ANY);
      using Owner = CODI_DD(T_Owner, CODI_ANY);
      using Impl = CODI_DD(T_Impl, CODI_IMPLEMENTATION);

      CODI_INLINE StaticThreadLocalPointerInterface() {}
      ~StaticThreadLocalPointerInterface() {}

      static CODI_INLINE void set(Type* other);
      static CODI_INLINE Type* get() const;
  };
}
