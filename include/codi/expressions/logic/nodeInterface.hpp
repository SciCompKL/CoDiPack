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

#include "../../config.h"
#include "../../misc/macros.hpp"

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
   * @tparam T_Impl
   */
  template<typename T_Impl>
  struct NodeInterface {
    public:

      using Impl = CODI_DD(T_Impl, NodeInterface);  ///< See NodeInterface.

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
