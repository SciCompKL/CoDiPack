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
   * @brief Computes the binomial coefficent n over k.
   *
   * @param[in] n  The set size n.
   * @param[in] k  The selection size k.
   *
   * @return THe binomial coefficent n over k.
   */
  CODI_INLINE size_t binomial(size_t n, size_t k) {
    if(k == 0) {
      return 1;
    } else if(n < k) {
      return 0; // outsize of the domain we assume zero values.
    } else if(n == k) {
      return 1;
    } else {
      return binomial(n - 1, k - 1) + binomial(n - 1, k);
    }
  }

  /**
   * @brief A helper namespace that computes the binomial coefficent of n over k at compile time.
   */
  namespace BinomalTemplate {

    /* Forward declaration of the binomial class */
    template<int n, int k>
    struct Binomial;

    /**
     * @brief Selector for the binomial if the requested data is outside of the definition.
     *
     * This class is specialized for the parameter outOfBounds, which selects then the
     * normal computation or the the out of bounds handling.
     *
     * @tparam           n  The set size n.
     * @tparam           k  The selection size k.
     * @tparam outOfBounds  The indicator if n < l
     */
    template<int n, int k, bool outOfBounds>
    struct BinomialCompute;

    /**
     * @brief Specialization for the binomal for the out of bounds case.
     *
     * The value is defined as zero.
     *
     * @tparam           n  The set size n.
     * @tparam           k  The selection size k.
     */
    template<int n, int k>
    struct BinomialCompute<n, k, true> {

      /** @brief The out of bounds value is zero */
      static const int value = 0;
    };

    /**
     * @brief Specialization for the binomal for the regular computation.
     *
     * The value is computed via the addition formulation of binomail coefficent.
     *
     * @tparam           n  The set size n.
     * @tparam           k  The selection size k.
     */
    template<int n, int k>
    struct BinomialCompute<n, k, false> {

      /** @brief The addition formulation of the binomial coefficent. */
      static const int value = (Binomial<n-1,k-1>::value + Binomial<n-1,k>::value);
    };

    /**
     * @brief Binomal value definition and computation during compile time.
     *
     * The value is computed by the helper class BinomalCompute, that selects if the
     * provided arguments are out of bounds or not.
     *
     * The value is computed via the addition formulation of binomail coefficent.
     *
     * @tparam           n  The set size n.
     * @tparam           k  The selection size k.
     */
    template<int n, int k>
    struct Binomial {
      /** @brief The addition formulation of the binomial coefficent. */
      static const int value =  BinomialCompute<n, k, n < k>::value;
    };

    /**
     * @brief Specialization for the 0, 0 case.
     *
     * This is necessary to avoid an ambigious error.
     */
    template<>
    struct Binomial<0,0> {
      /** @brief The zero zero case defines a value of one */
      static const int value = 1;
    };

    /**
     * @brief Specialization for the n, 0 case.
     *
     * @tparam           n  The set size n.
     */
    template<int n>
    struct Binomial<n,0> {
      /** @brief The n zero case defines a value of one */
      static const int value = 1;
    };

    /**
     * @brief Specialization for the n, n case.
     *
     * @tparam           n  The set size n.
     */
    template<int n>
    struct Binomial<n,n> {
      /** @brief The n n case defines a value of one */
      static const int value = 1;
    };

    /**
     * @brief Computes the binomal coefficent n over k at compile time.
     *
     * @return The binomial coefficient as a compile time constant.
     *
     * @tparam n  The set size n.
     * @tparam k  The selection size k.
     */
    template<size_t n, size_t k>
    static const size_t binomial() {
      return Binomial<n,k>::value;
    }
  }

  /**
   * @brief Computes the binomal coefficent n over k at compile time.
   *
   * The binomial coefficient as a compile time constant is stored in value.
   *
   * @tparam n  The set size n.
   * @tparam k  The selection size k.
   */
  template<size_t n, size_t k>
  struct Binomial {
    /** @brief The binomial coefficient as a compile time constant. */
    static const size_t value = BinomalTemplate::Binomial<n,k>::value;
  };

  /**
   * @brief Computes the binomal coefficent n over k at compile time.
   *
   * @return The binomial coefficient as a compile time constant.
   *
   * @tparam n  The set size n.
   * @tparam k  The selection size k.
   */
  template<size_t n, size_t k>
  const size_t binomial() {
    return BinomalTemplate::Binomial<n,k>::value;
  }
}
