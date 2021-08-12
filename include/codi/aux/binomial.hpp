#pragma once

#include "../config.h"
#include "macros.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Binomial coefficient computation.
   *
   * Recursive implementation of
   *
   * \f[
   *  \binom{n}{k} = \frac{n!}{k!(n-k)!} = \binom{n-1}{k} + \binom{n-1}{k-1}
   * \f]
   *
   * @param[in] n  Number of choices.
   * @param[in] k  Number of draws.
   * @return Binomial coefficient.
   */
  CODI_INLINE size_t constexpr binomial(size_t n, size_t k) {
    // clang-format off
    return
        k == 0 ? 1 : (
        n < k  ? 0 : ( // Outside of the domain we assume zero values.
        n == k ? 1 : (
        /* default */ binomial(n - 1, k - 1) + binomial(n - 1, k)
        )));
    // clang-format on
  }
}
