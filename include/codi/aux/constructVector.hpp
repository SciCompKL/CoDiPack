#pragma once

#include <array>
#include <vector>

#include "../config.h"
#include "macros.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /// Vector construction helper
  ///
  /// @tparam A vector type
  template<typename _V, typename = void>
  struct ConstructVectorImpl {
    public:

      using V = CODI_DD(_V, CODI_T(std::vector<CODI_ANY>));  ///< See ConstructVectorImpl

      /// Default implementation calls one argument constructor.
      static V construct(size_t const size) {
        return V(size);
      }
  };

  /// Specialization for std::array
  ///
  /// @tparam _T  Any type
  /// @tparam _n  array size.
  template<typename _T, size_t _n>
  struct ConstructVectorImpl<std::array<_T, _n>> {
    public:

      using T = CODI_DD(_T, CODI_ANY);  ///< See ConstructVectorImpl
      static size_t constexpr n = _n;   ///< See ConstructVectorImpl

      /// Only asserts the argument for the correct size.
      static std::array<T, n> construct(size_t const size) {
        codiAssert(size == n);

        return std::array<T, n>();
      }
  };

  /// Helper for the construction of vector types provided by the user.
  template<typename V>
  V constructVector(size_t const size) {
    return ConstructVectorImpl<V>::construct(size);
  }
}
