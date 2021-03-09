#pragma once

#include <vector>

#include "../../aux/macros.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  template<typename _T>
  struct HessianInterface {
    public:

      using T = CODI_DECLARE_DEFAULT(_T, double);

      virtual ~HessianInterface() {}

      size_t getM() const;
      size_t getN() const;

      T operator()(size_t const i, size_t const j, size_t const k) const;
      T& operator()(size_t const i, size_t const j, size_t const k);

      void resize(size_t const m, size_t const n);
      size_t size() const;
  };
}
