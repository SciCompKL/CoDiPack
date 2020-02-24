#pragma once

#include "../aux/macros.h"
#include "../config.h"
#include "logic/nodeInterface.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  template<typename _Real, typename _Impl>
  struct ExpressionInterface {
    public:

      using Real = DECLARE_DEFAULT(_Real, double);
      using Impl = DECLARE_DEFAULT(_Impl, ExpressionInterface);

      using StoreAs = ExpressionInterface;

      CODI_INLINE const Impl& cast() const {
        return static_cast<const Impl&>(*this);
      }

      CODI_INLINE const Real getValue() const;

      template<size_t argNumber>
      CODI_INLINE Real getJacobian() const;

    private:
      ExpressionInterface& operator=(const ExpressionInterface&) = delete;
  };
}
