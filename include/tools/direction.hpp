/**
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015 Chair for Scientific Computing (SciComp), TU Kaiserslautern
 * Homepage: http://www.scicomp.uni-kl.de
 * Contact:  Prof. Nicolas R. Gauger (codi@scicomp.uni-kl.de)
 *
 * Lead developers: Max Sagebaum, Tim Albring (SciComp, TU Kaiserslautern)
 *
 * This file is part of CoDiPack (http://www.scicomp.uni-kl.de/software/codi).
 *
 * CoDiPack is free software: you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation, either version 2 of the
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

#include "../configure.h"

/**
 * @brief Global namespace for CoDiPack - Code Differentiation Package
 */
namespace codi {

  /**
   * @brief The vector for direction of the forward mode or the reverse mode.
   *
   * The direction is an array with a static size that can be added and multiplied by
   * a scalar value.
   *
   * @tparam Real  The scalar value type that is used by the array.
   * @tparam  dim  The dimension of the array.
   */
  template<typename Real, size_t dim>
  class Direction {
    private:
      Real vector[dim]; /**< The data vector with the given dimension */

    public:

      /**
       * @brief Creates a zero direction.
       */
      CODI_INLINE Direction() :
        vector() {}

      /**
       * @brief Get the i-th element of the direction.
       *
       * No bounds checks are performed.
       *
       * @param[in] i  The index for the vector.
       *
       * @return The element from the vector.
       */
      CODI_INLINE Real& operator[] (const size_t& i) {
        return vector[i];
      }

      /**
       * @brief Get the i-th element of the direction.
       *
       * No bounds checks are performed.
       *
       * @param[in] i  The index for the vector.
       *
       * @return The element from the vector.
       */
      CODI_INLINE const Real& operator[] (const size_t& i) const {
        return vector[i];
      }

      /**
       * @brief Assign operator for the direction.
       *
       * @param[in] v  The values from the direction are set to the values of this direction object.
       *
       * @return Reference to this object.
       */
      CODI_INLINE Direction<Real, dim>& operator = (const Direction<Real, dim>& v) {
        for(size_t i = 0; i < dim; ++i) {
          this->vector[i] = v.vector[i];
        }

        return *this;
      }

      /**
       * @brief Update operator for the direction.
       *
       * @param[in] v  The values from the direction are added to the values of this direction object.
       *
       * @return Reference to this object.
       */
      CODI_INLINE Direction<Real, dim>& operator += (const Direction<Real, dim>& v) {
        for(size_t i = 0; i < dim; ++i) {
          this->vector[i] += v.vector[i];
        }

        return *this;
      }
  };

  /**
   * @brief Scalar multiplication of a direction.
   *
   * Performs the operation w = s * v
   *
   * @param[in] s  The scalar value for the multiplication.
   * @param[in] v  The direction that is multiplied.
   *
   * @return The direction with the result.
   *
   * @tparam Real  The scalar value type that is used by the direction.
   * @tparam  dim  The dimension of the direction.
   */
  template<typename Real, size_t dim>
  CODI_INLINE Direction<Real, dim> operator * (const Real& s, const Direction<Real, dim>& v) {
    Direction<Real, dim> r;
    for(size_t i = 0; i < dim; ++i) {
      r[i] = s * v[i];
    }

    return r;
  }

  /**
   * @brief Scalar multiplication of a direction.
   *
   * Performs the operation w = v * s
   *
   * @param[in] v  The direction that is multiplied.
   * @param[in] s  The scalar value for the multiplication.
   *
   * @return The direction with the result.
   *
   * @tparam Real  The scalar value type that is used by the direction.
   * @tparam  dim  The dimension of the direction.
   */
  template<typename Real, size_t dim>
  CODI_INLINE Direction<Real, dim> operator * (const Direction<Real, dim>& v, const Real& s) {
    return s * v;
  }

  /**
   * @brief Check if at least one component of the direction is not equal to s.
   *
   * The operator returns false if all the components of the direction v are equal to s,
   * true otherwise.
   *
   * @param[in] s  The scalar value that is checked against the components of the direction
   * @param[in] v  The direction that is compared with the scalar value.
   *
   * @return true if at least one component of v is not equal to s, false otherwise.
   *
   * @tparam    A  The type of the scalr value.
   * @tparam Real  The scalar value type that is used by the direction.
   * @tparam  dim  The dimension of the direction.
   */
  template<typename A, typename Real, size_t dim>
  CODI_INLINE bool operator != (const A& s, const Direction<Real, dim>& v) {
    for(size_t i = 0; i < dim; ++i) {
      if( s != v[i] ) {
        return true;
      }
    }

    return false;
  }

  /**
   * @brief Check if at least one component of the direction is not equal to s.
   *
   * The operator returns false if all the components of the direction v are equal to s,
   * true otherwise.
   *
   * @param[in] v  The direction that is compared with the scalar value.
   * @param[in] s  The scalar value that is checked against the components of the direction
   *
   * @return true if at least one component of v is not equal to s, false otherwise.
   *
   * @tparam    A  The type of the scalr value.
   * @tparam Real  The scalar value type that is used by the direction.
   * @tparam  dim  The dimension of the direction.
   */
  template<typename A, typename Real, size_t dim>
  CODI_INLINE bool operator != (const Direction<Real, dim>& v, const A& s) {
    return s != v;
  }
}
