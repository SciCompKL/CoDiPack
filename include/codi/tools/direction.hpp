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

#include <initializer_list>

#include "../configure.h"
#include "../typeFunctions.hpp"

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
       * @brief Creates a direction with the same value in every component.
       *
       * @param[in] s  The value that is set to all components.
       */
      CODI_INLINE Direction(const Real& s) :
        vector()
      {
        for(size_t i = 0; i < dim; ++i) {
          vector[i] = s;
        }
      }

      /**
       * @brief The direction is initialized with the values from the initializer list.
       *
       * If the list is to small, then only the first m elements are set.
       * If the list is to large, then only the first dim elements are set.
       *
       * @param[in] l  The list with the values for the direction.
       */
      CODI_INLINE Direction(std::initializer_list<Real> l) :
        vector()
      {
        size_t size = std::min(dim, l.size());
        const Real* array = l.begin(); // this is possible because the standard requires an array storage
        for(size_t i = 0; i < size; ++i) {
          vector[i] = array[i];
        }
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

      /**
       * @brief Checks if all entries in the direction are also a total zero.
       *
       * @return true if all entries are a total zero.
       */
      CODI_INLINE bool isTotalZero() const {
        for(size_t i = 0; i < dim; ++i) {
          if( !codi::isTotalZero(vector[i])) {
            return false;
          }
        }

        return true;
      }
  };

  /**
   * @brief Tests if all elements of the given direction are finite.
   *
   * Calls on all elements codi::isfinite.
   *
   * @param[in] d  The direction vector that is tested.
   * @return true if all elements are finite.
   * @tparam Real  The computation type of the direction vector.
   * @tparam  dim  The dimension of the direction vector.
   */
  template<typename Real, size_t dim>
  bool isfinite(const Direction<Real, dim>& d) {
    bool finite = true;

    for(size_t i = 0; i < dim; ++i) {
      finite &= codi::isfinite(d[i]);
    }

    return finite;
  }

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
   * @brief Scalar division of a direction.
   *
   * Performs the operation w = v / s
   *
   * @param[in] v  The direction that is divided.
   * @param[in] s  The scalar value for the division.
   *
   * @return The direction with the result.
   *
   * @tparam Real  The scalar value type that is used by the direction.
   * @tparam  dim  The dimension of the direction.
   */
  template<typename Real, size_t dim>
  CODI_INLINE Direction<Real, dim> operator / (const Direction<Real, dim>& v, const Real& s) {
    Direction<Real, dim> r;
    for(size_t i = 0; i < dim; ++i) {
      r[i] = v[i] / s;
    }

    return r;
  }

  /**
   * @brief Addition of two directions.
   *
   * Performs the operation w = v1 + v2
   *
   * @param[in] v1  The first direction that is added.
   * @param[in] v2  The second direction that is added.
   *
   * @return The direction with the result.
   *
   * @tparam Real  The scalar value type that is used by the direction.
   * @tparam  dim  The dimension of the direction.
   */
  template<typename Real, size_t dim>
  CODI_INLINE Direction<Real, dim> operator + (const Direction<Real, dim>& v1, const Direction<Real, dim>& v2) {
    Direction<Real, dim> r;
    for(size_t i = 0; i < dim; ++i) {
      r[i] = v1[i] + v2[i];
    }

    return r;
  }

  /**
   * @brief Subtraction of two directions.
   *
   * Performs the operation w = v1 - v2
   *
   * @param[in] v1  The first direction that is added.
   * @param[in] v2  The second direction that is subtracted.
   *
   * @return The direction with the result.
   *
   * @tparam Real  The scalar value type that is used by the direction.
   * @tparam  dim  The dimension of the direction.
   */
  template<typename Real, size_t dim>
  CODI_INLINE Direction<Real, dim> operator - (const Direction<Real, dim>& v1, const Direction<Real, dim>& v2) {
    Direction<Real, dim> r;
    for(size_t i = 0; i < dim; ++i) {
      r[i] = v1[i] - v2[i];
    }

    return r;
  }

  /**
   * @brief Negation of a direction.
   *
   * Performs the negation on all elements.
   *
   * @param[in] v  The first direction that is added.
   *
   * @return The direction with the result.
   *
   * @tparam Real  The scalar value type that is used by the direction.
   * @tparam  dim  The dimension of the direction.
   */
  template<typename Real, size_t dim>
  CODI_INLINE Direction<Real, dim> operator - (const Direction<Real, dim>& v) {
    Direction<Real, dim> r;
    for(size_t i = 0; i < dim; ++i) {
      r[i] = -v[i];
    }

    return r;
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
   * @tparam    A  The type of the scalar value.
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
   * @tparam    A  The type of the scalar value.
   * @tparam Real  The scalar value type that is used by the direction.
   * @tparam  dim  The dimension of the direction.
   */
  template<typename A, typename Real, size_t dim>
  CODI_INLINE bool operator != (const Direction<Real, dim>& v, const A& s) {
    return s != v;
  }

  /**
   * @brief Output the direction to a stream.
   *
   * The output format is: {v[0], v[1], ..., v[dim - 1]}
   *
   * @param[in,out] os  The output stream that is used for the writing.
   * @param[in]      v  The direction that is written to the stream.
   *
   * @return The output stream os.
   *
   * @tparam Real  The scalar value type that is used by the direction.
   * @tparam  dim  The dimension of the direction.
   */
  template<typename Real, size_t dim>
  std::ostream& operator<<(std::ostream& os, const Direction<Real, dim>& v){
    os << "{";
    for(size_t i = 0; i < dim; ++i) {
      if(i != 0) {
        os << ", ";
      }
      os << v[i];
    }
    os << "}";

    return os;
  }
}
