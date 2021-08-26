#pragma once

#include "../../aux/macros.hpp"
#include "../../config.h"

/** \copydoc codi::Namespace */
namespace codi {

  /// Helper class for the delayed write access of a reference.
  ///
  /// This class can be returned instead of a reference when the owner of the reference wants to be informed about write
  /// actions to the reference. Each assign call is forwarded to `data.setLogic(i,j, v)`.
  ///
  /// @tparam T_Impl  The issuing class of the delay accessor.
  template<typename T_Impl>
  struct JacobianDelayAccessor {
    public:

      using Impl = CODI_DD(T_Impl, CODI_ANY);  ///< See DelayAccessor.

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
