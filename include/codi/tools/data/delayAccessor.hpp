#pragma once

#include "../../aux/macros.hpp"
#include "../../config.h"

/** \copydoc codi::Namespace */
namespace codi {
  template<typename _Impl>
  struct DelayAccessor {
    public:

      using Impl = CODI_DD(_Impl, CODI_ANY);

    private:

      size_t i;
      size_t j;

      Impl& data;

    public:

      DelayAccessor(size_t const i, size_t const j, Impl& data) : i(i), j(j), data(data) {}

      template<typename T>
      DelayAccessor& operator=(T const& v) {
        data.setLogic(i, j, v);

        return *this;
      }

      operator typename Impl::T() const {
        return const_cast<Impl const&>(data).operator()(i, j);
      }
  };
}
