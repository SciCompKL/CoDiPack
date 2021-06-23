#pragma once

#include "../../aux/macros.hpp"
#include "../../config.h"

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

      using Impl = CODI_DD(_Impl, NodeInterface);  ///< See NodeInterface.

      /// Cast to the implementation.
      CODI_INLINE Impl const& cast() const {
        return static_cast<Impl const&>(*this);
      }

      /*******************************************************************************/
      /// @name Interface definition

      static bool constexpr EndPoint = CODI_UNDEFINED_VALUE;  ///< If this expression is handled as a leaf in the tree.

      /// Call the link method of the given logic for all arguments (links) of this node (not to be confused with args).
      /// Pass args to each call.
      template<typename Logic, typename... Args>
      CODI_INLINE void forEachLink(TraversalLogic<Logic>& logic, Args&&... args) const;

      /// Call the link method of the given logic for all arguments (links) of this node (not to be confused with args).
      /// Pass args to each call.
      template<typename Logic, typename... Args>
      CODI_INLINE static typename Logic::ResultType constexpr forEachLinkConstExpr(Args&&... args);
  };
}
