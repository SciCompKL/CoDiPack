/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2021 Chair for Scientific Computing (SciComp), TU Kaiserslautern
 * Homepage: http://www.scicomp.uni-kl.de
 * Contact:  Prof. Nicolas R. Gauger (codi@scicomp.uni-kl.de)
 *
 * Lead developers: Max Sagebaum, Tim Albring (SciComp, TU Kaiserslautern)
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
 * Authors:
 *  - SciComp, TU Kaiserslautern:
 *     Max Sagebaum
 *     Tim Albring
 *     Johannes Blühdorn
 */

#pragma once

#include "configure.h"
#include "macros.h"
#include "tools/direction.hpp"

/**
 * @brief Global namespace for CoDiPack - Code Differentiation Package
 */
namespace codi {

  /**
   * @brief Gradient value traits are used to access gradient data in a general way.
   *
   * The functions in these traits provide a generalized access to gradient values
   * which can be used to write code that can handle arbitrary gradient values.
   *
   * @tparam T An arbitrary data type. Should be a floating point type.
   */
  template<typename T>
  struct GradientValueTraits {
      using Data = T; /**< Inner data type of the gradient value. Default: Same as T */

      /**
       * @brief Get the vector size of the data.
       *
       * @return Default: 1
       */
      static CODI_INLINE constexpr size_t getVectorSize() {
        return 1;
      }

      /**
       * @brief Access an entry of the gradient value as read write reference.
       *
       * @param[in,out]   v  The gradient value
       * @param[in]     pos  The array access position.
       *                     Range: [0 ,..., getVectorSize())
       *
       * @return The element at the specified position.
       */
      static CODI_INLINE T& at(T& v, const size_t pos) {
        CODI_UNUSED(pos);

        return v;
      }

      /**
       * @brief Access an entry of the gradient value as read reference.
       *
       * \copydetails at()
       */
      static CODI_INLINE const T& at(const T& v, const size_t pos) {
        CODI_UNUSED(pos);

        return v;
      }
  };

  /**
   * @brief Specialization of GradientValueTraits for Direction.
   */
  template<typename T, size_t n>
  struct GradientValueTraits<Direction<T, n>> {

      using Data = T; /**< Entry type of the direction */

      /**
       * @brief The array size of the Direction.
       *
       * @return n
       */
      static CODI_INLINE constexpr size_t getVectorSize() {
        return n;
      }

      /**
       * \copydoc GradientValueTraits::at()
       */
      static CODI_INLINE T& at(Direction<T,n>& v, const size_t pos) {
        return v[pos];
      }

      /**
       * \copydoc GradientValueTraits::at(const T& v, const size_t pos)
       */
      static CODI_INLINE const T& at(const Direction<T,n>& v, const size_t pos) {
        return v[pos];
      }
  };
}
