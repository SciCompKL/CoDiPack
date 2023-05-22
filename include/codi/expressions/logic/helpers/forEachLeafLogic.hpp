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

#include <type_traits>
#include <utility>

#include "../../../config.h"
#include "../../../misc/macros.hpp"
#include "../../../traits/expressionTraits.hpp"
#include "../traversalLogic.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Implement logic for leaf nodes only.
   *
   * This class calls:
   *  - handleActive for each leaf node that implements LhsExpressionInterface,
   *  - handleConstant for each leaf node that implements ConstantExpression.
   *
   * For details about the expression traversal see TraversalLogic.
   *
   * @tparam T_Impl  Class implementing this interface.
   */
  template<typename T_Impl>
  struct ForEachLeafLogic : public TraversalLogic<T_Impl> {
    public:

      using Impl = CODI_DD(T_Impl, CODI_T(TraversalLogic<CODI_ANY>));  ///< See ForEachLeafLogic.

      /*******************************************************************************/
      /// @name Interface definition
      /// @{

      /// Called for leaf nodes which implement LhsExpressionInterface.
      template<typename Node, typename... Args>
      void handleActive(Node const& node, Args&&... args) {
        CODI_UNUSED(node, args...);
      }

      /// Called for leaf nodes which implement ConstantExpression.
      template<typename Node, typename... Args>
      void handleConstant(Node const& node, Args&&... args) {
        CODI_UNUSED(node, args...);
      }

      /// @}
      /*******************************************************************************/
      /// @name Overwrites from TraversalLogic
      /// @{

      /// \copydoc codi::TraversalLogic::leaf()
      template<typename Node, typename... Args>
      CODI_INLINE ExpressionTraits::EnableIfLhsExpression<Node> leaf(Node const& node, Args&&... args) {
        cast().handleActive(node, std::forward<Args>(args)...);
      }

      /// \copydoc codi::TraversalLogic::leaf()
      template<typename Node, typename... Args>
      CODI_INLINE ExpressionTraits::EnableIfConstantExpression<Node> leaf(Node const& node, Args&&... args) {
        cast().handleConstant(node, std::forward<Args>(args)...);
      }

      /// @}

    private:

      CODI_INLINE Impl& cast() {
        return static_cast<Impl&>(*this);
      }
  };
}
