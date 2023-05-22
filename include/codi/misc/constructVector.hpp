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

#include <array>
#include <vector>

#include "../config.h"
#include "macros.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /// Vector construction helper.
  ///
  /// @tparam A vector type
  template<typename T_V, typename = void>
  struct ConstructVectorImpl {
    public:

      using V = CODI_DD(T_V, CODI_T(std::vector<CODI_ANY>));  ///< See ConstructVectorImpl.

      /// Default implementation assumes that there is a constructor that takes the vector size as its single argument.
      static V construct(size_t const size) {
        return V(size);
      }
  };

  /// Specialization for std::array.
  ///
  /// @tparam T_T  Any type.
  /// @tparam T_n  Array size.
  template<typename T_T, size_t T_n>
  struct ConstructVectorImpl<std::array<T_T, T_n>> {
    public:

      using T = CODI_DD(T_T, CODI_ANY);  ///< See ConstructVectorImpl
      static size_t constexpr n = T_n;   ///< See ConstructVectorImpl

      /// Only asserts the argument for the correct size.
      static std::array<T, n> construct(size_t const size) {
        CODI_UNUSED(size);
        codiAssert(size == n);

        return std::array<T, n>();
      }
  };

  /// Helper for the construction of vector types provided by the user.
  template<typename V>
  V constructVector(size_t const size) {
    return ConstructVectorImpl<V>::construct(size);
  }
}
