/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2017 Chair for Scientific Computing (SciComp), TU Kaiserslautern
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
 * Authors: Max Sagebaum, Tim Albring, (SciComp, TU Kaiserslautern)
 */

#pragma once

#include <type_traits>

/**
 * @brief Global namespace for CoDiPack - Code Differentiation Package
 */
namespace codi {

#ifndef DOXYGEN_DISABLE
  // default definition
  template <typename T, const bool isArithmetic>
  struct IsTotalZeroImpl {};

  //call t == 0 for all arithmetic types e.g. double, int @internal */
  template <typename T>
  struct IsTotalZeroImpl<T, true> {
      static CODI_INLINE bool isTotalZero(const T &t) {
        return t == T();
      }
  };

  //call t.isTotalZero for all other types */
  template <typename T>
  struct IsTotalZeroImpl<T, false> {
      static CODI_INLINE bool isTotalZero(const T &t) {
        return t.isTotalZero();
      }
  };
#endif

  /**
   * @brief Check if all values in the type are zero
   *
   * The function is used to determine if a Jacobian of that type should be stored on the tape.
   * It is also used to determine if the adjoint update should be performed.
   *
   * On arithmetic types the implementation calls t == 0.0.
   * On other types the implementation calls t.isTotalZero();
   *
   * @param[in] t  The value that is checked for a total zero.
   * @return true if all primal values and gradient values(if they exist are zero) otherwise false.
   */
  template <typename T>
  CODI_INLINE bool isTotalZero(const T& t) {
    return IsTotalZeroImpl<T, std::is_arithmetic<T>::value>::isTotalZero(t);}
}
