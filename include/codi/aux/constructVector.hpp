#pragma once

#include <array>
#include <vector>

#include "../config.h"
#include "macros.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /// Vector construction helper.
  ///
  /// @tparam A vector type
  template<typename T_V, typename = void>
  struct ConstructVectorImpl {
    public:

      using V = CODI_DD(T_V, CODI_T(std::vector<CODI_ANY>));  ///< See ConstructVectorImpl.

      /// Default implementation assumes that there is a constructor that takes the vector size as its single argument.
      static V construct(size_t const size) {
        return V(size);
      }
  };

  /// Specialization for std::array.
  ///
  /// @tparam T_T  Any type.
  /// @tparam T_n  Array size.
  template<typename T_T, size_t T_n>
  struct ConstructVectorImpl<std::array<T_T, T_n>> {
    public:

      using T = CODI_DD(T_T, CODI_ANY);  ///< See ConstructVectorImpl
      static size_t constexpr n = T_n;   ///< See ConstructVectorImpl

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
