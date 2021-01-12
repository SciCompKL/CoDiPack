#pragma once

/** \copydoc codi::Namespace */
namespace codi {

  struct DummyValue {
    public:
      template<typename T>
      CODI_INLINE void operator=(T const& v) {
        CODI_UNUSED(v);
      }
  };

  struct DummyJacobian {
    public:
      CODI_INLINE DummyValue operator()(size_t const i, size_t const j) const {
        CODI_UNUSED(i);
        CODI_UNUSED(j);

        return DummyValue();
      }
  };
}
