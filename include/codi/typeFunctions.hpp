/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2019 Chair for Scientific Computing (SciComp), TU Kaiserslautern
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
#include "expressions.hpp"

/**
 * @brief Global namespace for CoDiPack - Code Differentiation Package
 */
namespace codi {

#ifndef DOXYGEN_DISABLE
  // default definition. call t.isTotalZero()
  template <typename T, typename Enable = void>
  struct IsTotalZeroImpl {
      static CODI_INLINE bool isTotalZero(const T &t) {
        return t.isTotalZero();
      }
  };

  //call t == 0 for all arithmetic types e.g. double, int @internal
  template <typename T>
  struct IsTotalZeroImpl<
      T,
      typename std::enable_if<std::is_arithmetic<T>::value>::type
      >
  {
      static CODI_INLINE bool isTotalZero(const T &t) {
        return t == T();
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
    return IsTotalZeroImpl<T>::isTotalZero(t);}

#ifndef DOXYGEN_DISABLE
  // Take address of a T instance
  template <typename T, typename Enable = void>
  struct AddressOfImpl {
      typedef typename std::add_pointer<T>::type PointerType;

      static CODI_INLINE PointerType get(T &t) {
          return &t;
      }
  };
#endif

  /**
   * @brief Return address of a variable
   *
   * The default implementation returns &t.
   *
   * @param[in] t The value from which the adress is taken.
   * @tparam T Type of the variable.
   */
  template <typename T>
  CODI_INLINE
  typename AddressOfImpl<T>::PointerType addressof(T& t) {
    return AddressOfImpl<T>::get(t);
  }

#ifndef DOXYGEN_DISABLE
  // check if variable is finite
  template <typename T, typename Enable = void>
  struct IsFiniteImpl {
      static CODI_INLINE bool get(const T &t) {
          using std::isfinite;
          return isfinite(t);
      }
  };

  //call the specialized isfinite implementation for codi::Expression
  template <typename T>
  struct IsFiniteImpl<
    T,
    typename std::enable_if<
      std::is_base_of<
        codi::Expression<typename codi::TypeTraits<T>::Real, T>,
        T
      >::value
    >::type
  >
  {
      static CODI_INLINE bool get(const T &t) {
        using std::isfinite;
        return isfinite(t.getValue());
      }
  };

#endif

  /**
   * @brief Check if variable is finite
   *
   * The default implementation calls isfinite(t) without a namespace specifier.
   *
   * @param[in] t The value for which the is finite attribute is checked.
   * @tparam T Type of the variable
   */
  template <typename T>
  CODI_INLINE bool isfinite(const T& t) {
    return IsFiniteImpl<T>::get(t);
  }

#ifndef DOXYGEN_DISABLE
  // check if variable is finite
  template <typename T>
  struct ArrayAccessImpl {
      static CODI_INLINE size_t get(const T &t) {
          return t;
      }
  };
#endif

  /**
   * @brief Helper function to convert a CoDiPack index into a size_t object.
   *
   * The default implementation converts T to size_t.
   *
   * @param[in] t The value which is converted.
   * @return The converted value.
   *
   * @tparam T Type of the variable
   */
  template <typename T>
  CODI_INLINE size_t arrayAccess(const T& t) {
    return ArrayAccessImpl<T>::get(t);
  }
}
