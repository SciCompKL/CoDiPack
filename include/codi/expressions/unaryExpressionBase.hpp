/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2022 Chair for Scientific Computing (SciComp), TU Kaiserslautern
 * Homepage: http://www.scicomp.uni-kl.de
 * Contact:  Prof. Nicolas R. Gauger (codi@scicomp.uni-kl.de)
 *
 * Lead developers: Max Sagebaum, Johannes Blühdorn (SciComp, TU Kaiserslautern)
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
 *  - SciComp, TU Kaiserslautern:
 *    - Max Sagebaum
 *    - Johannes Blühdorn
 *    - Former members:
 *      - Tim Albring
 */
#pragma once

#include "../config.h"
#include "../misc/macros.hpp"
#include "expressionInterface.hpp"
#include "logic/compileTimeTraversalLogic.hpp"
#include "logic/nodeInterface.hpp"
#include "logic/traversalLogic.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Base implementation for expressions with on argument.
   *
   * Implements the NodeInterface and provides the member for the expression argument. The functions getValue and
   * getJacobian from the ExpressionInterface have to be implemented by the derived class.
   *
   * See \ref Expressions "Expression" design documentation for details about the expression system in CoDiPack.
   *
   * @tparam T_Real  Original primal value of the statement/expression.
   * @tparam T_Arg  The ExpressionInterface type of the argument.
   * @tparam T_Impl  Class implementing this interface.
   */
  template<typename T_Real, typename T_Arg, typename T_Impl>
  struct UnaryExpressionBase : public ExpressionInterface<T_Real, T_Impl> {
    public:

      using Real = CODI_DD(T_Real, double);                                       ///< See UnaryExpressionBase.
      using Arg = CODI_DD(T_Arg, CODI_T(ExpressionInterface<double, CODI_ANY>));  ///< See UnaryExpressionBase.
      using Impl = CODI_DD(T_Impl, ExpressionInterface);                          ///< See UnaryExpressionBase.

      typename Arg::StoreAs arg;  ///< Argument of the expression.

      /// Constructor
      template<typename RealArg>
      CODI_INLINE explicit UnaryExpressionBase(ExpressionInterface<RealArg, Arg> const& arg) : arg(arg.cast()) {}

      CODI_INLINE UnaryExpressionBase(UnaryExpressionBase const&) = default;  ///< Constructor

      /*******************************************************************************/
      /// @name Implementation of ExpressionInterface
      /// @{

      using StoreAs = Impl;                   ///< \copydoc codi::ExpressionInterface::StoreAs
      using ADLogic = typename Arg::ADLogic;  ///< \copydoc codi::ExpressionInterface::ADLogic

      /// @}
      /*******************************************************************************/
      /// @name Implementation of NodeInterface
      /// @{

      static bool constexpr EndPoint = false;  ///< \copydoc codi::NodeInterface::EndPoint

      /// \copydoc codi::NodeInterface::forEachLink
      template<typename Logic, typename... Args>
      CODI_INLINE void forEachLink(TraversalLogic<Logic>& logic, Args&&... args) const {
        logic.cast().template link<0>(arg, this->cast(), std::forward<Args>(args)...);
      }

      /// \copydoc codi::NodeInterface::forEachLinkConstExpr
      template<typename CompileTimeLogic, typename... Args>
      CODI_INLINE static typename CompileTimeLogic::ResultType constexpr forEachLinkConstExpr(Args&&... args) {
        return CompileTimeLogic::template link<0, Arg, Impl>(std::forward<Args>(args)...);
      }

      /// @}
  };
}
