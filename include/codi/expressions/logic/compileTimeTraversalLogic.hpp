#pragma once

#include <utility>

#include "../../config.h"
#include "../../aux/macros.hpp"
#include "nodeInterface.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  template<typename _ResultType, typename _Impl>
  struct CompileTimeTraversalLogic {
    public:

      using ResultType = CODI_DECLARE_DEFAULT(_ResultType, size_t);
      using Impl = CODI_DECLARE_DEFAULT(_Impl, CompileTimeTraversalLogic);

      static constexpr ResultType NeutralElement = {};

      CODI_INLINE static constexpr ResultType reduce(ResultType a, ResultType b) {
        return a + b;
      }

      template<typename Node, typename ... Args>
      CODI_INLINE static constexpr ResultType node(Args&& ... args) {
        // Default logic forwards to all links
        return toLinks<Node>(std::forward<Args>(args)...);
      }

      template<typename Node, typename ... Args>
      CODI_INLINE static constexpr ResultType term(Args&& ... CODI_UNUSED_ARG(args)) {
        // Default logic does nothing
        return Impl::NeutralElement;
      }

      template<size_t LeafNumber, typename Leaf, typename Root, typename ... Args>
      CODI_INLINE static constexpr ResultType link(Args&& ... args) {
        // Default logic forwards to node evaluation
        return toNode<Leaf>(std::forward<Args>(args)...);
      }

      template<typename Node, typename ... Args>
      CODI_INLINE static constexpr ResultType eval(Args&& ... args) {
        return toNode<Node>(std::forward<Args>(args)...);
      }

    protected:

      template<typename TraversalImpl, bool endPoint = false>
      struct CallSwitch {
          template<typename Node, typename ... Args>
          CODI_INLINE static constexpr ResultType call(Args&& ... args) {
            return TraversalImpl::template node<Node>(std::forward<Args>(args)...);
          }
      };

      template<typename TraversalImpl>
      struct CallSwitch <TraversalImpl, true>{
          template<typename Node, typename ... Args>
          CODI_INLINE static constexpr ResultType call(Args&& ... args) {
            return TraversalImpl::template term<Node>(std::forward<Args>(args)...);
          }
      };

      template<typename Node, typename ... Args>
      CODI_INLINE static constexpr ResultType toNode(Args&& ... args) {
        return CallSwitch<Impl, Node::EndPoint>::template call<Node>(std::forward<Args>(args)...);
      }

      template<typename Node, typename ... Args>
      CODI_INLINE static constexpr ResultType toLinks(Args&& ... args) {
        return Node::template forEachLinkConst<Impl>(std::forward<Args>(args)...);
      }
  };
}
