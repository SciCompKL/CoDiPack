#pragma once

#include <utility>

#include "../../aux/macros.hpp"
#include "../../config.h"
#include "nodeInterface.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Traversal of CoDiPack expressions during compile time.
   *
   * For a detailed explanation of the traversal structure please see \ref customExpressionLogic "Expression traversal".
   *
   * All information needs to be provided as arguments and all computations must be constexpr.
   *
   * @tparam _ResultType  Type of the computed result.
   * @tparam _Impl  Class implementing this interface.
   */
  template<typename _ResultType, typename _Impl>
  struct CompileTimeTraversalLogic {
    public:

      using ResultType = CODI_DD(_ResultType, size_t);         ///< See CompileTimeTraversalLogic
      using Impl = CODI_DD(_Impl, CompileTimeTraversalLogic);  ///< See CompileTimeTraversalLogic

      static constexpr ResultType NeutralElement = {};  ///< Neutral element of the reduction.

      /**
       * @brief Start the evaluation of the logic on the given expression.
       */
      template<typename Node, typename... Args>
      CODI_INLINE static constexpr ResultType eval(Args&&... args) {
        return toNode<Node>(std::forward<Args>(args)...);
      }

      /*******************************************************************************/
      /// @name Interface definition
      /// @{

      /**
       * @brief Reduction operation for the results of two links.
       *
       * Needs to be a constexpr.
       *
       * Default: summation
       */
      CODI_INLINE static constexpr ResultType reduce(ResultType a, ResultType b) {
        return a + b;
      }

      /**
       * @brief Called for each node in the expression.
       *
       * Implementations can call the toLinks method in order to evaluate all links of the node.
       *
       * Needs to be a constexpr.
       *
       * Default: Call each link of the node and forward all arguments.
       */
      template<typename Node, typename... Args>
      CODI_INLINE static constexpr ResultType node(Args&&... args) {
        // Default logic forwards to all links
        return toLinks<Node>(std::forward<Args>(args)...);
      }

      /**
       * @brief Called for all termination nodes in the expression.
       *
       * Needs to be a constexpr.
       *
       * Default: Returns NeutralElement
       */
      template<typename Node, typename... Args>
      CODI_INLINE static constexpr ResultType term(Args&&... CODI_UNUSED_ARG(args)) {
        // Default logic does nothing
        return Impl::NeutralElement;
      }

      /**
       * @brief Called for all links in the expression.
       *
       * Implementations can call the toNode method in order to evaluate either term or node depending on the leaf.
       *
       * Needs to be a constexpr.
       *
       * Default: Call the leaf node and forward all arguments.
       */
      template<size_t LeafNumber, typename Leaf, typename Root, typename... Args>
      CODI_INLINE static constexpr ResultType link(Args&&... args) {
        // Default logic forwards to node evaluation
        return toNode<Leaf>(std::forward<Args>(args)...);
      }

    /// @}

    protected:

    #ifndef DOXYGEN_DISABLE
      template<typename TraversalImpl, bool endPoint = false>
      struct CallSwitch {
        template<typename Node, typename... Args>
        CODI_INLINE static constexpr ResultType call(Args&&... args) {
          return TraversalImpl::template node<Node>(std::forward<Args>(args)...);
        }
      };

      template<typename TraversalImpl>
      struct CallSwitch<TraversalImpl, true> {
        template<typename Node, typename... Args>
        CODI_INLINE static constexpr ResultType call(Args&&... args) {
          return TraversalImpl::template term<Node>(std::forward<Args>(args)...);
        }
      };
      #endif

      /// Helper method to distinguish between termination nodes and normal nodes.
      template<typename Node, typename... Args>
      CODI_INLINE static constexpr ResultType toNode(Args&&... args) {
        return CallSwitch<Impl, Node::EndPoint>::template call<Node>(std::forward<Args>(args)...);
      }

      /// Helper method which calls forEachLinkConstExpr on the node.
      template<typename Node, typename... Args>
      CODI_INLINE static constexpr ResultType toLinks(Args&&... args) {
        return Node::template forEachLinkConstExpr<Impl>(std::forward<Args>(args)...);
      }
  };
}
