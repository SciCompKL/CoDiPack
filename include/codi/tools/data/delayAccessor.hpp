#pragma once

#include "../../aux/macros.hpp"
#include "../../config.h"

/** \copydoc codi::Namespace */
namespace codi {

  /// Helper class for the delayed access of a reference.
  ///
  /// The class will forward assign calls to `data.setLogic(i,j, v)`.
  ///
  /// @tparam _Impl  The issuing class of the delay accessor.
  template<typename _Impl>
  struct JacobianDelayAccessor {
    public:

      using Impl = CODI_DD(_Impl, CODI_ANY);  ///< See DelayAccessor.

    private:

      size_t i;
      size_t j;

      Impl& data;

    public:

      /// Constructor
      JacobianDelayAccessor(size_t const i, size_t const j, Impl& data) : i(i), j(j), data(data) {}

      /// Forwards to `data.setLogic(i, j, v)`.
      template<typename T>
      JacobianDelayAccessor& operator=(T const& v) {
        data.setLogic(i, j, v);

        return *this;
      }

      /// Convert to the underlying type.
      operator typename Impl::T() const {
        return const_cast<Impl const&>(data).operator()(i, j);
      }
  };
}
