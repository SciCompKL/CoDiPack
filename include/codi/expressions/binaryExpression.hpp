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
#include "../traits/expressionTraits.hpp"
#include "expressionInterface.hpp"
#include "logic/compileTimeTraversalLogic.hpp"
#include "logic/nodeInterface.hpp"
#include "logic/traversalLogic.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Interface for implementing the logic for a BinaryExpression.
   *
   * See \ref Expressions "Expression" design documentation for details about the expression system in CoDiPack.
   *
   * @tparam T_Real  Original primal value of the statement/expression.
   */
  template<typename T_Real>
  struct BinaryOperation {
    public:

      using Real = CODI_DD(T_Real, double);  ///< See BinaryOperation.

      /// Compute the primal value from the arguments.
      ///
      /// The type of the arguments is the result of a getValue call on the expressions.
      template<typename ArgA, typename ArgB>
      static CODI_INLINE Real primal(ArgA const& argA, ArgB const& argB);

      /// Compute the gradient with respect to the first argument
      ///
      /// The type of the arguments is the type of the result of a getValue call on the expressions.
      template<typename ArgA, typename ArgB>
      static CODI_INLINE Real gradientA(ArgA const& argA, ArgB const& argB, Real const& result);

      /// Compute the gradient with respect to the second argument
      ///
      /// The type of the arguments is the type of the result of a getValue call on the expressions.
      template<typename ArgA, typename ArgB>
      static CODI_INLINE Real gradientB(ArgA const& argA, ArgB const& argB, Real const& result);
  };

  /**
   * @brief Represents an operator with two arguments in the expression tree.
   *
   * See \ref Expressions "Expression" design documentation for details about the expression system in CoDiPack.
   *
   * @tparam T_Real  Original primal value of the statement/expression.
   * @tparam T_ArgA  The ExpressionInterface type of the first argument.
   * @tparam T_ArgB  The ExpressionInterface type of the second argument.
   * @tparam T_Operation  The logic for computing the primal value and Jacobians. Must implement BinaryOperation.
   */
  template<typename T_Real, typename T_ArgA, typename T_ArgB, template<typename> class T_Operation>
  struct BinaryExpression : public ExpressionInterface<T_Real, BinaryExpression<T_Real, T_ArgA, T_ArgB, T_Operation> > {
    public:
      using Real = CODI_DD(T_Real, double);                                                 ///< See BinaryExpression.
      using ArgA = CODI_DD(T_ArgA, CODI_T(ExpressionInterface<double, CODI_ANY>));          ///< See BinaryExpression.
      using ArgB = CODI_DD(T_ArgB, CODI_T(ExpressionInterface<double, CODI_ANY>));          ///< See BinaryExpression.
      using Operation = CODI_DD(CODI_T(T_Operation<Real>), CODI_T(BinaryOperation<Real>));  ///< See BinaryExpression.

      typename ArgA::StoreAs argA;  ///< First argument of the expression.
      typename ArgB::StoreAs argB;  ///< Second argument of the expression.
      Real result;                  ///< Precomputed result.

      /// Constructor
      template<typename RealA, typename RealB>
      explicit BinaryExpression(ExpressionInterface<RealA, ArgA> const& argA,
                                ExpressionInterface<RealB, ArgB> const& argB)
          : argA(argA.cast()),
            argB(argB.cast()),
            result(Operation::primal(this->argA.getValue(), this->argB.getValue())) {}

      /*******************************************************************************/
      /// @name Implementation of ExpressionInterface
      /// @{

      using StoreAs = BinaryExpression;  ///< \copydoc codi::ExpressionInterface::StoreAs
      using ActiveResult =
          typename ExpressionTraits::ValidateResult<typename ArgA::ActiveResult, typename ArgB::ActiveResult>::
              ActiveResult;  ///< \copydoc codi::ExpressionInterface::ActiveResult

      /// \copydoc codi::ExpressionInterface::getValue()
      CODI_INLINE Real const& getValue() const {
        return result;
      }

      /// \copydoc codi::ExpressionInterface::getJacobian()
      template<size_t argNumber>
      CODI_INLINE Real getJacobian() const {
        if (0 == argNumber) {
          return Operation::gradientA(argA.getValue(), argB.getValue(), result);
        } else {
          return Operation::gradientB(argA.getValue(), argB.getValue(), result);
        }
      }

      /// @}
      /*******************************************************************************/
      /// @name Implementation of NodeInterface
      /// @{

      static bool constexpr EndPoint = false;  ///< \copydoc codi::NodeInterface::EndPoint

      /// \copydoc codi::NodeInterface::forEachLink
      template<typename Logic, typename... Args>
      CODI_INLINE void forEachLink(TraversalLogic<Logic>& logic, Args&&... args) const {
        logic.cast().template link<0>(argA, *this, std::forward<Args>(args)...);
        logic.cast().template link<1>(argB, *this, std::forward<Args>(args)...);
      }

      /// \copydoc codi::NodeInterface::forEachLinkConstExpr
      template<typename CompileTimeLogic, typename... Args>
      CODI_INLINE static typename CompileTimeLogic::ResultType constexpr forEachLinkConstExpr(Args&&... args) {
        return CompileTimeLogic::reduce(
            CompileTimeLogic::template link<0, ArgA, BinaryExpression>(std::forward<Args>(args)...),
            CompileTimeLogic::template link<1, ArgB, BinaryExpression>(std::forward<Args>(args)...));
      }

      /// @}
  };
}
