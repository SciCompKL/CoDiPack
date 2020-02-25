#pragma once

#include "../config.h"
#include "../aux/macros.h"

/** \copydoc codi::Namespace */
namespace codi {

  struct TraversalLogic;

  template<typename _Impl>
  struct NodeInterface {
    public:

      using Impl = DECLARE_DEFAULT(_Impl, NodeInterface);

      static bool constexpr EndPoint;

      CODI_INLINE Impl const& cast() const {
        return static_cast<Impl const&>(*this);
      }

      template<typename Logic, typename ... Args>
      CODI_INLINE void forEachLink(TraversalLogic<Logic>& logic, Args&& ... args) const;

      template<typename Logic, typename ... Args>
      CODI_INLINE static typename Logic::ResultType constexpr forEachLinkConst(Args&& ... args);
  };
}
