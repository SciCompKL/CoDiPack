#pragma once

#include <utility>

#include "../../config.h"
#include "../../aux/macros.hpp"
#include "nodeInterface.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  template<typename _Impl>
  struct TraversalLogic {
    public:

      using Impl = CODI_DECLARE_DEFAULT(_Impl, TraversalLogic);

      CODI_INLINE Impl& cast() {
        return static_cast<Impl&>(*this);
      }

      template<typename Node, typename ... Args>
      CODI_INLINE void node(Node const& node, Args&& ... args) {
        // Default logic forwards to all links
        toLinks(node, std::forward<Args>(args)...);
      }

      template<typename Node, typename ... Args>
      CODI_INLINE void term(Node const& node, Args&& ... args) {
        CODI_UNUSED(node, args...);
        // Default logic does nothing
      }

      template<size_t LeafNumber, typename Leaf, typename Root, typename ... Args>
      CODI_INLINE void link(Leaf const& leaf, Root const& root, Args&& ... args) {
        CODI_UNUSED(root, args...);
        // Default logic forwards to node evaluation
        toNode(leaf, std::forward<Args>(args)...);
      }

      template<typename Node, typename ... Args>
      CODI_INLINE void eval(NodeInterface<Node> const& node, Args&& ... args) {
        toNode(node.cast(), std::forward<Args>(args)...);
      }

    protected:

      template<typename TraversalImpl, bool endPoint = false>
      struct CallSwitch {
          template<typename ... Args>
          CODI_INLINE static void call(TraversalImpl& impl, Args&& ... args) {
            impl.node(std::forward<Args>(args)...);
          }
      };

      template<typename TraversalImpl>
      struct CallSwitch <TraversalImpl, true>{
          template<typename ... Args>
          CODI_INLINE static void call(TraversalImpl& impl, Args&& ... args) {
            impl.term(std::forward<Args>(args)...);
          }
      };

      template<typename Node, typename ... Args>
      CODI_INLINE void toNode(Node const& node, Args&& ... args) {
        CallSwitch<Impl, Node::EndPoint>::call(cast(), node, std::forward<Args>(args)...);
      }

      template<typename Node, typename ... Args>
      CODI_INLINE void toLinks(Node const& node, Args&& ... args) {
        node.forEachLink(cast(), std::forward<Args>(args)...);
      }
  };
}

