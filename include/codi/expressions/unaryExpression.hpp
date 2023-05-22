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

#include "../config.h"
#include "../misc/macros.hpp"
#include "expressionInterface.hpp"
#include "logic/compileTimeTraversalLogic.hpp"
#include "logic/nodeInterface.hpp"
#include "logic/traversalLogic.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Interface for implementing the logic for a UnaryExpression.
   *
   * See \ref Expressions "Expression" design documentation for details about the expression system in CoDiPack.
   *
   * @tparam T_Real  Original primal value of the statement/expression.
   */
  template<typename T_Real>
  struct UnaryOperation {
    public:

      using Real = CODI_DD(T_Real, double);  ///< See UnaryOperation.

      /// Compute the primal value from the argument.
      ///
      /// The type of the argument is the type of the result of a getValue call on the expression.
      template<typename Arg>
      static CODI_INLINE Real primal(Arg const& arg);

      /// Compute the gradient with respect to the argument
      ///
      /// The type of the argument is the type of the result of a getValue call on the expression.
      template<typename Arg>
      static CODI_INLINE Real gradient(Arg const& arg, Real const& result);
  };

  /**
   * @brief Represents an operator with one argument in the expression tree.
   *
   * See \ref Expressions "Expression" design documentation for details about the expression system in CoDiPack.
   *
   * @tparam T_Real  Original primal value of the statement/expression.
   * @tparam T_Arg  The ExpressionInterface type of the argument.
   * @tparam T_Operation  The logic for computing the primal value and Jacobian. Must implement UnaryOperation.
   */
  template<typename T_Real, typename T_Arg, template<typename> class T_Operation>
  struct UnaryExpression : public ExpressionInterface<T_Real, UnaryExpression<T_Real, T_Arg, T_Operation> > {
    public:

      using Real = CODI_DD(T_Real, double);                                                ///< See UnaryExpression.
      using Arg = CODI_DD(T_Arg, CODI_T(ExpressionInterface<double, CODI_ANY>));           ///< See UnaryExpression.
      using Operation = CODI_DD(CODI_T(T_Operation<Real>), CODI_T(UnaryOperation<Real>));  ///< See UnaryExpression.

      using ActiveResult = typename Arg::ActiveResult;  ///< See ExpressionInterface.

      typename Arg::StoreAs arg;  ///< Argument of the expression.
      Real result;                ///< Precomputed result.

      /// Constructor
      template<typename RealArg>
      explicit UnaryExpression(ExpressionInterface<RealArg, Arg> const& arg)
          : arg(arg.cast()), result(Operation::primal(this->arg.getValue())) {}

      /*******************************************************************************/
      /// @name Implementation of ExpressionInterface
      /// @{

      using StoreAs = UnaryExpression;  ///< \copydoc codi::ExpressionInterface::StoreAs

      /// \copydoc codi::ExpressionInterface::getValue()
      CODI_INLINE Real const& getValue() const {
        return result;
      }

      /// \copydoc codi::ExpressionInterface::getJacobian()
      template<size_t argNumber>
      CODI_INLINE Real getJacobian() const {
        return Operation::gradient(arg.getValue(), result);
      }

      /// @}
      /*******************************************************************************/
      /// @name Implementation of NodeInterface
      /// @{

      static bool constexpr EndPoint = false;  ///< \copydoc codi::NodeInterface::EndPoint

      /// \copydoc codi::NodeInterface::forEachLink
      template<typename Logic, typename... Args>
      CODI_INLINE void forEachLink(TraversalLogic<Logic>& logic, Args&&... args) const {
        logic.cast().template link<0>(arg, *this, std::forward<Args>(args)...);
      }

      /// \copydoc codi::NodeInterface::forEachLinkConstExpr
      template<typename CompileTimeLogic, typename... Args>
      CODI_INLINE static typename CompileTimeLogic::ResultType constexpr forEachLinkConstExpr(Args&&... args) {
        return CompileTimeLogic::template link<0, Arg, UnaryExpression>(std::forward<Args>(args)...);
      }

      /// @}
  };
}
