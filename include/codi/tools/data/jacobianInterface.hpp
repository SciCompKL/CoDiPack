#pragma once

#include "../../aux/macros.hpp"
#include "../../config.h"

/** \copydoc codi::Namespace */
namespace codi {

template<typename _T>
  struct JacobianInterface {

      using T = CODI_DECLARE_DEFAULT(_T, double);

      size_t getM() const;
      size_t getN() const;

      T operator()(size_t const i, size_t const j) const;

      T& operator()(size_t const i, size_t const j);

      void resize(size_t const m, size_t const n);
      size_t size() const;
  };
}
