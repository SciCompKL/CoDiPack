#pragma once

#include "../../config.h"
#include "../../aux/macros.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  template<typename Logic>
  struct TraversalLogic;

  /**
   * @brief Node side interface for the traversal of expressions.
   *
   * See the TraversalLogic and CompileTimeTraversalLogic for details on how this interface is used.
   *
   * Implementations need to call the link methods for each argument of the node.
   *
   * @tparam _Impl
   */
  template<typename _Impl>
  struct NodeInterface {
    public:

      using Impl = CODI_DD(_Impl, NodeInterface);  ///< See NodeInterface

      /// Cast to the implementation.
      CODI_INLINE Impl const& cast() const {
        return static_cast<Impl const&>(*this);
      }

      /*******************************************************************************/
      /// @name Interface definition

      static bool constexpr EndPoint = CODI_UNDEFINED_VALUE;  ///< If this expression is handled as a termination point.

      /// Call the link method on the logic for all arguments of this node.
      template<typename Logic, typename... Args>
      CODI_INLINE void forEachLink(TraversalLogic<Logic>& logic, Args&&... args) const;

      /// Call the link method on the logic for all arguments of this node.
      template<typename Logic, typename... Args>
      CODI_INLINE static typename Logic::ResultType constexpr forEachLinkConstExpr(Args&&... args);
  };
}
