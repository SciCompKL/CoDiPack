#pragma once

#include <utility>

#include "../../aux/macros.hpp"
#include "../../config.h"
#include "nodeInterface.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Traversal of CoDiPack expressions.
   *
   * For a detailed explanation of the traversal structure please see \ref customExpressionLogic "Expression traversal".
   *
   * Implementing classes can have members which can store information required for the traversal.
   *
   * @tparam _Impl  Class implementing this interface.
   */
  template<typename _Impl>
  struct TraversalLogic {
    public:

      using Impl = CODI_DD(_Impl, TraversalLogic);  ///< See TraversalLogic

      /// Cast to the implementation.
      CODI_INLINE Impl& cast() {
        return static_cast<Impl&>(*this);
      }

      /**
       * @brief Start the evaluation of the logic on the given expression.
       */
      template<typename Node, typename... Args>
      CODI_INLINE void eval(NodeInterface<Node> const& node, Args&&... args) {
        toNode(node.cast(), std::forward<Args>(args)...);
      }

      /*******************************************************************************/
      /// @name Interface definition
      /// @{

      /**
       * @brief Called for each node in the expression.
       *
       * Implementations can call the toLinks method in order to evaluate all links of the node.
       *
       * Default: Call each link of the node and forward all arguments.
       */
      template<typename Node, typename... Args>
      CODI_INLINE void node(Node const& node, Args&&... args) {
        // Default logic forwards to all links
        toLinks(node, std::forward<Args>(args)...);
      }

      /**
       * @brief Called for all termination nodes in the expression.
       *
       * Default: Does nothing.
       */
      template<typename Node, typename... Args>
      CODI_INLINE void term(Node const& node, Args&&... args) {
        CODI_UNUSED(node, args...);
        // Default logic does nothing
      }

      /**
       * @brief Called for all links in the expression.
       *
       * Implementations can call the toNode method in order to evaluate either term or node depending on the leaf.
       *
       * Default: Call the leaf node and forward all arguments.
       */
      template<size_t LeafNumber, typename Leaf, typename Root, typename... Args>
      CODI_INLINE void link(Leaf const& leaf, Root const& root, Args&&... args) {
        CODI_UNUSED(root, args...);
        // Default logic forwards to node evaluation
        toNode(leaf, std::forward<Args>(args)...);
      }

      /// @}

    protected:

#ifndef DOXYGEN_DISABLE
      template<typename TraversalImpl, bool endPoint = false>
      struct CallSwitch {
        public:
          template<typename... Args>
          CODI_INLINE static void call(TraversalImpl& impl, Args&&... args) {
            impl.node(std::forward<Args>(args)...);
          }
      };

      template<typename TraversalImpl>
      struct CallSwitch<TraversalImpl, true> {
        public:
          template<typename... Args>
          CODI_INLINE static void call(TraversalImpl& impl, Args&&... args) {
            impl.term(std::forward<Args>(args)...);
          }
      };
#endif

      /// Helper method to distinguish between termination nodes and normal nodes.
      template<typename Node, typename... Args>
      CODI_INLINE void toNode(Node const& node, Args&&... args) {
        CallSwitch<Impl, Node::EndPoint>::call(cast(), node, std::forward<Args>(args)...);
      }

      /// Helper method which calls forEachLink on the node.
      template<typename Node, typename... Args>
      CODI_INLINE void toLinks(Node const& node, Args&&... args) {
        node.forEachLink(cast(), std::forward<Args>(args)...);
      }
  };
}
