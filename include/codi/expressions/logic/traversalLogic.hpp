#pragma once

#include "../config.h"
#include "../aux/macros.h"
#include "nodeInterface.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  template<typename _Impl>
  struct TraversalLogic {
    public:

      using Impl = DECLARE_DEFAULT(_Impl, TraversalLogic);

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
        // Default logic does nothing
      }

      template<typename Leaf, typename Root, size_t LeafNumber, typename ... Args>
      CODI_INLINE void link(Leaf const& leaf, Root const& root, Args&& ... args) {
        // Default logic forwards to all links
        toNode(leaf, std::forward<Args>(args)...);
      }

      template<typename Node, typename ... Args>
      CODI_INLINE void eval(NodeInterface<Real, Node> const& node, Args&& ... args) {
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

      template<typename Impl, typename ... Args>
      CODI_INLINE void toNode(NodeInterface<Impl> const& node, Args&& ... args) {
        CallSwitch<Impl, Impl::EndPoint>::call(cast(), node, std::forward<Args>(args)...);
      }

      template<typename Node, typename ... Args>
      CODI_INLINE void toLinks(Node const& node, Args&& ... args) {
        node.forEachLink(cast(), std::forward<Args>(args)...);
      }
  };
}

