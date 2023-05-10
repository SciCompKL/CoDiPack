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
   * @brief Traversal of CoDiPack expressions.
   *
   * For a detailed explanation of the traversal structure please see \ref customExpressionLogic "Expression traversal".
   *
   * Implementing classes can have members which store information required for the traversal.
   *
   * @tparam T_Impl  Class implementing this interface.
   */
  template<typename T_Impl>
  struct TraversalLogic {
    public:

      using Impl = CODI_DD(T_Impl, TraversalLogic);  ///< See TraversalLogic.

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
        // Default logic forwards to all links.
        toLinks(node, std::forward<Args>(args)...);
      }

      /**
       * @brief Called for all leaf nodes in the expression.
       *
       * Default: Does nothing.
       */
      template<typename Node, typename... Args>
      CODI_INLINE void leaf(Node const& node, Args&&... args) {
        CODI_UNUSED(node, args...);
        // Default logic does nothing.
      }

      /**
       * @brief Called for all links in the expression.
       *
       * Implementations can call the toNode method in order to evaluate either leaf or node depending on the child.
       *
       * Default: Call the child node and forward all arguments.
       */
      template<size_t ChildNumber, typename Child, typename Root, typename... Args>
      CODI_INLINE void link(Child const& child, Root const& root, Args&&... args) {
        CODI_UNUSED(root, args...);
        // Default logic forwards to node evaluation.
        toNode(child, std::forward<Args>(args)...);
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
            impl.leaf(std::forward<Args>(args)...);
          }
      };
#endif

      /// Helper method to distinguish between leaf nodes and normal nodes.
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
