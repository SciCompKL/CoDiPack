/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2025 Chair for Scientific Computing (SciComp), University of Kaiserslautern-Landau
 * Homepage: http://scicomp.rptu.de
 * Contact:  Prof. Nicolas R. Gauger (codi@scicomp.uni-kl.de)
 *
 * Lead developers: Max Sagebaum, Johannes Blühdorn (SciComp, University of Kaiserslautern-Landau)
 *
 * This file is part of CoDiPack (http://scicomp.rptu.de/software/codi).
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

      static size_t constexpr LinkCount = CODI_UNDEFINED_VALUE;  ///< Number of links the expression has. Zero will
                                                                 ///< handle the expression as a leaf node.

      template<size_t argNumber>
      CODI_INLINE CODI_UNDEFINED const& getLink() const;  ///< Get the argument for the specific link. Usually
                                                                  ///< an expression.
  };
}
