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

#include <vector>
#include <array>

#include "../../configure.h"
#include "../../macros.h"

/**
 * @brief Global namespace for CoDiPack - Code Differentiation Package
 */
namespace codi {

  /**
   *  @brief Adapters for arrays such that they can be provided as template template arguments.
   */
  namespace adapters {

    /**
     * @brief Template template argument adapters for std::array
     *
     * @tparam s  The size of the array.
     */
    template<size_t s>
    struct StdArray {

        /**
         * @brief Type definition for std::array that can be used as a template template argument.
         *
         * @tparam T  The storage type of the array.
         */
        template <typename T>
        using Type = std::array<T, s>;
    };

    /**
     * @brief Type definition for std::vector that can be used as a template template argument.
     *
     * @tparam T  The storage type of the vector.
     */
    template <typename T>
    using StdVector = std::vector<T>;

  }

  /**
   * @brief Empty definition of a vector storage object
   *
   * This interface captures the access to an arbitrary data store. In CoDiPack it is specialized
   * for std::vector and std::array. If a user wants to use the Jacobian or Hessian data structures with other
   * data objects in the background this class needs to be specialized.
   *
   * @tparam Vec  The implementation of the vector storage.
   */
  template <typename Vec>
  struct VectorStorage {
      using VecType = Vec; /**< The type this specialization is implemented for. */
      using Element = char; /**< The type of the data entries for the storage. */

      /**
       * @brief Constructor that initializes the storage with a given size.
       *
       * @param[in] size  The size of the data.
       */
      explicit VectorStorage(size_t size);

      /**
       * @brief Get access to the data of the storage.
       *
       * @return A pointer to the data representation.
       */
      CODI_INLINE Element const* data() const;

      /** \copydoc data() const */
      CODI_INLINE Element* data();

      /**
       * @brief Reference to the i-th data element.
       * @param[in] i  The position in the data vector.
       * @return The reference to the i-th element.
       */
      CODI_INLINE Element& operator[](size_t i);

      /**
       * @brief The current size of the storage data.
       * @return The size of the storage data.
       */
      CODI_INLINE size_t size();

      /**
       * @brief Resize the underlying data to fit the new size.
       *
       * If the vector type does not allow to be resized, then an exception is thrown.
       * @param[in] size  The new size.
       */
      CODI_NO_INLINE void resize(size_t const size);
  };

  /**
   * @brief Specialization of the VectorStorage interface for std::vector
   *
   * All calls are forwarded to the equivalent routines in the std::vector interface.
   *
   * @tparam T  The elements in the std::vector
   * @tparam Allocator  The allocation object for the std::vector
   */
  template <typename T, typename Allocator>
  struct VectorStorage<std::vector<T, Allocator>> {

      using VecType = std::vector<T, Allocator>; /**< std::vector */
      using Element = T; /**< Element type of the std::vector */

      VecType vec; /**< Instantiation of the std::vector */

      /** \copydoc VectorStorage::VectorStorage */
      explicit VectorStorage(size_t size) : vec(size) {}

      /** \copydoc VectorStorage::data */
      CODI_INLINE T const* data() const {
        return vec.data();
      }

      /** \copydoc VectorStorage::data */
      CODI_INLINE T* data() {
        return vec.data();
      }

      /** \copydoc VectorStorage::operator[] */
      CODI_INLINE T& operator[](size_t i) {
        return vec[i];
      }

      /** \copydoc VectorStorage::size */
      CODI_INLINE size_t size() {
        return vec.size();
      }

      /**
       * \copybrief
       *
       * @param[in] size  The new size of the vector.
       */
      CODI_NO_INLINE void resize(size_t const size) {
        vec.resize(size);
      }
  };

  /**
   * @brief Specialization of the VectorStorage interface for std::array
   *
   * The resize function in this specialization throws an error if the new size is different then the template
   * parameter.
   *
   * All other calls are forwarded to the equivalent routines in the std::array interface.
   *
   * @tparam T  The elements in the std::array
   * @tparam N  The size of the std::array.
   */
  template <typename T, size_t N>
  struct VectorStorage<std::array<T, N>> {

      using VecType = std::array<T, N>; /**< std::array */
      using Element = T; /**< Element type of std::array */

      VecType vec; /**< Instantiation of the std::array */

      /**
       * @brief The constructor ignores the given size.
       *
       * @param[in] size  Unused
       */
      explicit VectorStorage(size_t size) : vec() {CODI_UNUSED(size);}

      /** \copydoc VectorStorage::data */
      CODI_INLINE T const* data() const {
        return vec.data();
      }

      /** \copydoc VectorStorage::data */
      CODI_INLINE T* data() {
        return vec.data();
      }

      /** \copydoc VectorStorage::operator[] */
      CODI_INLINE T& operator[](size_t i) {
        return vec[i];
      }

      /** \copydoc VectorStorage::size */
      CODI_INLINE size_t size() {
        return N;
      }

      /**
       * @brief Throws an exception if the given size is different than N
       *
       * @param[in] size  Unused
       */
      CODI_INLINE void resize(size_t const size) {
        if(N != size) {
          CODI_EXCEPTION("Can not resize std::array.");
        }
      }
  };
}
