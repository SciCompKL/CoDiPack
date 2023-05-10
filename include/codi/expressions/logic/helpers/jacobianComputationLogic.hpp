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
   *
   * @tparam T_Impl  Class implementing this interface.
   */
  template<typename T_Impl>
  struct JacobianComputationLogic : public TraversalLogic<T_Impl> {
    public:

      using Impl = CODI_DD(T_Impl, CODI_T(TraversalLogic<CODI_ANY>));  ///< See JacobianComputationLogic.

      /*******************************************************************************/
      /// @name Interface definition
      /// @{

      /// Called for leaf nodes which implement LhsExpressionInterface.
      template<typename Node, typename Jacobian, typename... Args>
      void handleJacobianOnActive(Node const& node, Jacobian jacobian, Args&&... args);

      /// @}
      /*******************************************************************************/
      /// @name Overwrites from TraversalLogic
      /// @{

      /// \copydoc codi::TraversalLogic::leaf()
      template<typename Node, typename Jacobian, typename... Args>
      CODI_INLINE ExpressionTraits::EnableIfLhsExpression<Node> leaf(Node const& node, Jacobian jacobian,
                                                                     Args&&... args) {
        cast().handleJacobianOnActive(node, jacobian, std::forward<Args>(args)...);
      }

      using TraversalLogic<Impl>::leaf;

      /// Computes the \ref sec_reverseAD "reverse" AD equation for this link.
      ///
      /// The Jacobian is multiplied with the Jacobian of the link. The result is forwarded to the child.
      template<size_t ChildNumber, typename Child, typename Root, typename Jacobian, typename... Args>
      CODI_INLINE void link(Child const& child, Root const& root, Jacobian const& jacobian, Args&&... args) {
        auto linkJacobian = root.template getJacobian<ChildNumber>();
        auto curJacobian = linkJacobian * jacobian;

        cast().toNode(child, curJacobian, std::forward<Args>(args)...);
      }

      /// @}

    private:

      CODI_INLINE Impl& cast() {
        return static_cast<Impl&>(*this);
      }
  };
}
