#pragma once

#include "../../aux/macros.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  struct DummyValue {
    public:
      template<typename T>
      CODI_INLINE void operator=(T const& v) {
        CODI_UNUSED(v);
      }
  };

  struct DummyVector {
    public:
      CODI_INLINE DummyValue operator[](size_t const i) {
        CODI_UNUSED(i);
        return DummyValue();
      }

      CODI_INLINE size_t size() const {
        return (size_t)0;
      }
  };

  struct DummyHessian {
    public:
      CODI_INLINE DummyValue operator()(size_t const i, size_t const j, size_t const k) {
        CODI_UNUSED(i, j, k);

        return DummyValue();
      }
  };
}
