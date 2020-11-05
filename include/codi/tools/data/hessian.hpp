/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2020 Chair for Scientific Computing (SciComp), TU Kaiserslautern
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

#include "dummyValue.hpp"
#include "vectorStorage.hpp"
#include "../../configure.h"

/**
 * @brief Global namespace for CoDiPack - Code Differentiation Package
 */
namespace codi {

  /**
   * @brief Basic interface definition what the algorithms in the EvaluationHelper, Algorithms, etc. classes use to
   * access the data.
   *
   * @tparam T  The data type of the internal storage.
   */
  template<typename T>
  struct HessianInterface {

      virtual ~HessianInterface(){}

      /**
       * @brief Constant access to the specified element of the hessian.
       *
       * @param[in] i  Output value of the function. Range: [0, m)
       * @param[in] j  Input value of the function. First derivative direction. Range: [0,n) (Fastest running index)
       * @param[in] k  Input value of the function. Second derivative direction. Range: [0,n) (Slowest running index)
       *
       * @return Value/reference of the specified location.
       */
      virtual T operator()(const size_t i, const size_t j, const size_t k) const = 0;

      /**
       * @brief Access to the specified element of the hessian.
       *
       * \copydetails operator()(const size_t i, const size_t j, const size_t k) const
       */
      virtual T& operator()(const size_t i, const size_t j, const size_t k) = 0;

  };

  /**
   * @brief A dummy hessian that can be accessed via the call operator.
   */
  struct DummyHessian {

      /**
       * @brief Returns a dummy value.
       *
       * @param[in] i  Not used
       * @param[in] j  Not used
       * @param[in] k  Not used
       * @return A dummy value
       */
      CODI_INLINE DummyValue operator()(const size_t i, const size_t j, const size_t k) {
        CODI_UNUSED(i);
        CODI_UNUSED(j);
        CODI_UNUSED(k);

        return DummyValue();
      }
  };

  /**
   * @brief Default hessian implementation for algorithms in CoDiPack.
   *
   * The data layout of the hessian is described in the \ref FunctionDefinition "Algorithms" documentation.
   * All mathematical symbol names are described there.
   *
   * Structures that implement the same functions can be used a the same places where this hessian implementation
   * is used in CoDiPack.
   *
   * @tparam Vec  Either std::vector or std::array. For other types VectorStorage needs to be specialized.
   */
  template <typename Vec>
  struct Hessian {

    private:
      VectorStorage<Vec> values; /**< Storage for the hessian matrix */

      size_t m; /**< Number of function outputs */
      size_t n; /**< Number of function inputs */

    public:
      /**
       * @brief Inner element of the vector type.
       */
      using T = typename VectorStorage<Vec>::Element;

      /**
       * @brief Create a hessian with the given output and input size.
       *
       * @param[in] m  Number of function outputs
       * @param[in] n  Number of function inputs
       */
      explicit Hessian(size_t m, size_t n) : values(n * n * m), m(m), n(n) {}

      /**
       * @brief Get the number of function outputs.
       * @return Number of function outputs.
       */
      CODI_INLINE size_t getM() const {
        return m;
      }

      /**
       * @brief Get the number of function inputs.
       * @return Number of function inputs.
       */
      CODI_INLINE size_t getN() const {
        return n;
      }

      /**
       * @brief Constant access to the specified element of the hessian.
       *
       * @param[in] i  Output value of the function. Range: [0, m)
       * @param[in] j  Input value of the function. First derivative direction. Range: [0,n) (Fastest running index)
       * @param[in] k  Input value of the function. Second derivative direction. Range: [0,n) (Slowest running index)
       *
       * @return Value/reference of the specified location.
       */
      CODI_INLINE T operator()(const size_t i, const size_t j, const size_t k) const {
        return values.data()[computeIndex(i,j,k)];
      }

      /**
       * @brief Access to the specified element of the hessian.
       *
       * \copydetails operator()(const size_t i, const size_t j, const size_t k) const
       */
      CODI_INLINE T& operator()(const size_t i, const size_t j, const size_t k) {
        return values.data()[computeIndex(i,j,k)];
      }

      /**
       * @brief Set value in the hessian.
       *
       * @param[in] i  Output value of the function. Range: [0, m)
       * @param[in] j  Input value of the function. First derivative direction. Range: [0,n) (Fastest running index)
       * @param[in] k  Input value of the function. Second derivative direction. Range: [0,n) (Slowest running index)
       * @param[in] v  The value that is set into the hessian.
       */
      template<typename T>
      CODI_INLINE void set(const size_t i, const size_t j, const size_t k, const T& v) {
        values.data()[computeIndex(i,j,k)] = v;
      }

    private:

      /**
       * @brief Computes the offset into the data array.
       *
       * @param[in] i  n. Mid running index
       * @param[in] j  1. Fastest running index
       * @param[in] k  m*n. Slowest running index
       * @return
       */
      CODI_INLINE size_t computeIndex(const size_t i, const size_t j, const size_t k) const {
        return k * n * m + i * n + j;
      }
  };

}
