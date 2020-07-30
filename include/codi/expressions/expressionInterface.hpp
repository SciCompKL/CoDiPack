#pragma once

#include "../aux/macros.hpp"
#include "../config.h"
#include "logic/nodeInterface.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  template<typename _Real, typename _Impl>
  struct ExpressionInterface : public NodeInterface<_Impl> {
    public:

      using Real = CODI_DECLARE_DEFAULT(_Real, double);
      using Impl = CODI_DECLARE_DEFAULT(_Impl, ExpressionInterface);

      CODI_INLINE Impl const& cast() const {
        return static_cast<Impl const&>(*this);
      }

      /*******************************************************************************
       * Section: Start of interface definition
       *
       */

      using StoreAs = ExpressionInterface;

      CODI_INLINE Real const getValue() const;

      template<size_t argNumber>
      CODI_INLINE Real getJacobian() const;

    private:
      ExpressionInterface& operator=(ExpressionInterface const&) = delete;
  };
}
