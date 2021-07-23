#pragma once

#include "../../aux/macros.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /// Dummy value that can be assigned.
  struct DummyValue {
    public:

      /// Ignore the assignment.
      template<typename T>
      CODI_INLINE void operator=(T const& v) {
        CODI_UNUSED(v);
      }
  };

  /// Dummy vector that provides a dummy element access and size function.
  struct DummyVector {
    public:

      /// Return a dummy value.
      CODI_INLINE DummyValue operator[](size_t const i) {
        CODI_UNUSED(i);
        return DummyValue();
      }

      /// Vector is always zero size.
      CODI_INLINE size_t size() const {
        return (size_t)0;
      }
  };
}
