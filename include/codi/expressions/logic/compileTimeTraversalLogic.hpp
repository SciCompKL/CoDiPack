/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2023 Chair for Scientific Computing (SciComp), University of Kaiserslautern-Landau
 * Homepage: http://www.scicomp.uni-kl.de
 * Contact:  Prof. Nicolas R. Gauger (codi@scicomp.uni-kl.de)
 *
 * Lead developers: Max Sagebaum, Johannes Blühdorn (SciComp, University of Kaiserslautern-Landau)
 *
 * This file is part of CoDiPack (http://www.scicomp.uni-kl.de/software/codi).
 *
 * CoDiPack is free software: you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * CoDiPack is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty
 * of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 * See the GNU General Public License for more details.
 * You should have received a copy of the GNU
 * General Public License along with CoDiPack.
 * If not, see <http://www.gnu.org/licenses/>.
 *
 * For other licensing options please contact us.
 *
 * Authors:
 *  - SciComp, University of Kaiserslautern-Landau:
 *    - Max Sagebaum
 *    - Johannes Blühdorn
 *    - Former members:
 *      - Tim Albring
 */
#pragma once

#include <utility>

#include "../../config.h"
#include "../../misc/macros.hpp"
#include "nodeInterface.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Traversal of CoDiPack expressions during compile time.
   *
   * For a detailed explanation of the traversal structure please see \ref customExpressionLogic "Expression traversal".
   *
   * All information must be provided as arguments and all computations must be constexpr.
   *
   * @tparam T_ResultType  Type of the computed result.
   * @tparam T_Impl  Class implementing this interface.
   */
  template<typename T_ResultType, typename T_Impl>
  struct CompileTimeTraversalLogic {
    public:

      using ResultType = CODI_DD(T_ResultType, size_t);         ///< See CompileTimeTraversalLogic.
      using Impl = CODI_DD(T_Impl, CompileTimeTraversalLogic);  ///< See CompileTimeTraversalLogic.

      static ResultType constexpr NeutralElement = {};  ///< Neutral element of the reduction.

      /**
       * @brief Start the evaluation of the logic on the given expression.
       */
      template<typename Node, typename... Args>
      CODI_INLINE static ResultType constexpr eval(Args&&... args) {
        return toNode<Node>(std::forward<Args>(args)...);
      }

      /*******************************************************************************/
      /// @name Interface definition
      /// @{

      /**
       * @brief Reduction operation for the results of two links.
       *
       * Must be a constexpr.
       *
       * Default: summation.
       */
      CODI_INLINE static ResultType constexpr reduce(ResultType a, ResultType b) {
        return a + b;
      }

      /**
       * @brief Called for each node in the expression.
       *
       * Implementations can call the toLinks method in order to evaluate all links of the node.
       *
       * Must be a constexpr.
       *
       * Default: Call each link of the node and forward all arguments.
       */
      template<typename Node, typename... Args>
      CODI_INLINE static ResultType constexpr node(Args&&... args) {
        // Default logic forwards to all links.
        return toLinks<Node>(std::forward<Args>(args)...);
      }

      /**
       * @brief Called for all leaf nodes in the expression.
       *
       * Must be a constexpr.
       *
       * Default: Returns NeutralElement.
       */
      template<typename Node, typename... Args>
      CODI_INLINE static ResultType constexpr leaf(Args&&... CODI_UNUSED_ARG(args)) {
        // Default logic does nothing.
        return Impl::NeutralElement;
      }

      /**
       * @brief Called for all links in the expression.
       *
       * Implementations can call the toNode method in order to evaluate either leaf or node depending on the child.
       *
       * Must be a constexpr.
       *
       * Default: Call the child node and forward all arguments.
       */
      template<size_t ChildNumber, typename Child, typename Root, typename... Args>
      CODI_INLINE static ResultType constexpr link(Args&&... args) {
        // Default logic forwards to node evaluation.
        return toNode<Child>(std::forward<Args>(args)...);
      }

      /// @}

    protected:

#ifndef DOXYGEN_DISABLE
      template<typename TraversalImpl, bool endPoint = false>
      struct CallSwitch {
        public:
          template<typename Node, typename... Args>
          CODI_INLINE static ResultType constexpr call(Args&&... args) {
            return TraversalImpl::template node<Node>(std::forward<Args>(args)...);
          }
      };

      template<typename TraversalImpl>
      struct CallSwitch<TraversalImpl, true> {
        public:
          template<typename Node, typename... Args>
          CODI_INLINE static ResultType constexpr call(Args&&... args) {
            return TraversalImpl::template leaf<Node>(std::forward<Args>(args)...);
          }
      };
#endif

      /// Helper method to distinguish between leaf nodes and normal nodes.
      template<typename Node, typename... Args>
      CODI_INLINE static ResultType constexpr toNode(Args&&... args) {
        return CallSwitch<Impl, Node::EndPoint>::template call<Node>(std::forward<Args>(args)...);
      }

      /// Helper method which calls forEachLinkConstExpr on the node.
      template<typename Node, typename... Args>
      CODI_INLINE static ResultType constexpr toLinks(Args&&... args) {
        return Node::template forEachLinkConstExpr<Impl>(std::forward<Args>(args)...);
      }
  };
}
