#pragma once

#include "../config.h"

/** \copydoc codi::Namespace */
namespace codi {

  inline size_t constexpr binomial(size_t n, size_t k) {
    return
        k == 0 ? 1 : (
        n < k  ? 0 : ( // outside of the domain we assume zero values
        n == k ? 1 : (
        /* default */ binomial(n - 1, k - 1) + binomial(n - 1, k)
        )));
  }
}
