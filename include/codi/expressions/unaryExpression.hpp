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
#include "unaryExpressionBase.hpp"

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
      using Jacobian = Real;                 ///< Can be overwritten by the implementation.

      /// Compute the primal value from the argument.
      ///
      /// The type of the argument is the type of the result of a getValue call on the expression.
      template<typename Arg>
      static CODI_INLINE Real primal(Arg const& arg);

      /// Compute the gradient with respect to the argument
      ///
      /// The type of the argument is the type of the result of a getValue call on the expression.
      template<typename Arg>
      static CODI_INLINE Jacobian gradient(Arg const& arg, Real const& result);
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
  struct UnaryExpression : public UnaryExpressionBase<T_Real, T_Arg, UnaryExpression<T_Real, T_Arg, T_Operation> > {
    public:

      using Real = CODI_DD(T_Real, double);                                                ///< See UnaryExpression.
      using Arg = CODI_DD(T_Arg, CODI_T(ExpressionInterface<double, CODI_ANY>));           ///< See UnaryExpression.
      using Operation = CODI_DD(CODI_T(T_Operation<Real>), CODI_T(UnaryOperation<Real>));  ///< See UnaryExpression.

      using Base = UnaryExpressionBase<T_Real, T_Arg, UnaryExpression>;  ///< Abbreviation of base class.
      using Jacobian = typename Operation::Jacobian;                     ///< Jacobian defined by the operation.

      Real result;  ///< Precomputed result.

      /// Constructor
      template<typename RealArg>
      CODI_INLINE explicit UnaryExpression(ExpressionInterface<RealArg, Arg> const& arg)
          : Base(arg), result(Operation::primal(this->arg.getValue())) {}

      CODI_INLINE UnaryExpression(UnaryExpression const&) = default;  ///< Constructor

      /*******************************************************************************/
      /// @name Implementation of ExpressionInterface
      /// @{

      /// \copydoc codi::ExpressionInterface::getValue()
      CODI_INLINE Real const& getValue() const {
        return result;
      }

      /// \copydoc codi::ExpressionInterface::getJacobian()
      template<size_t argNumber>
      CODI_INLINE Jacobian getJacobian() const {
        return Operation::gradient(Base::arg.getValue(), result);
      }

      /// @}
  };
}
