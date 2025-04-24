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

#include <tuple>

#include "../config.h"
#include "../misc/macros.hpp"
#include "../traits/computationTraits.hpp"
#include "../traits/expressionTraits.hpp"
#include "expressionInterface.hpp"
#include "logic/compileTimeTraversalLogic.hpp"
#include "logic/nodeInterface.hpp"
#include "logic/traversalLogic.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Interface for implementing the logic for a ComputeExpression.
   *
   * Represents an operation \f$ w = f(x)\f$.
   *
   * See \ref Expressions "Expression" design documentation for details about the expression system in CoDiPack.
   *
   * @tparam T_Real  Original primal value of the statement/expression.
   */
  template<typename T_Real>
  struct ComputeOperation {
    public:

      using Real = CODI_DD(T_Real, double);  ///< See ComputeOperation.

      /// Compute the primal value from the arguments. E.g. evaluates \f$ w = f(x)\f$
      ///
      /// The argument types are the types of the result of a getValue call on the expression.
      template<typename... Args>
      static CODI_INLINE Real primal(Args const&... args);

      /// Apply the forward AD mode with respect to the given argument on the tangent and return the result.
      /// Computes \f$ \dot w = \frac{\d f}{\d x_{argNumber}} \dot x_{argNumber} \f$
      template<size_t argNumber, typename Tangent, typename... Args>
      static CODI_INLINE auto applyTangent(Tangent const& tangent, Real const& result, Args const&... args);

      /// Apply the reverse AD mode with respect to the given argument on the adjoint and return the result.
      /// Computes \f$ \bar x_{argNumber} = \frac{\d f}{\d x_{argNumber}}^T \bar w \f$
      template<size_t argNumber, typename Adjoint, typename... Args>
      static CODI_INLINE auto applyAdjoint(Adjoint const& adjoint, Real const& result, Args const&... args);

      /// Get the math symbol of the operation. E.g. `+` for operators and `pow()` for functions.
      static CODI_INLINE std::string getMathRep();
  };

  /**
   * @brief Implements ComputeOperation for one argument.
   *
   * Implementations need to define #primal, #applyTangentArg and #applyAdjointArg.
   *
   * See ComputeOperation for details.
   *
   * @tparam T_Real  Original primal value of the statement/expression.
   * @tparam T_Impl  Final implementation of the operation.
   */
  template<typename T_Real, typename T_Impl>
  struct UnaryOperation : public ComputeOperation<T_Real> {
    public:

      using Real = CODI_DD(T_Real, double);  ///< See UnaryOperation.

      /// \copydoc codi::ComputeOperation::primal
      template<typename Arg>
      static CODI_INLINE Real primal(Arg const& arg);

      /// Apply the forward AD mode with respect to the argument on the tangent and return the result.
      /// Computes \f$ \dot w = \frac{\d f}{\d x} \dot x \f$
      template<typename Tangent, typename Arg>
      static CODI_INLINE auto applyTangentArg(Tangent const& tangent, Real const& result, Arg const& arg);

      /// Apply the reverse AD mode with respect to the argument on the adjoint and return the result.
      /// Computes \f$ \bar x = \frac{\d f}{\d x}^T \bar w \f$
      template<typename Adjoint, typename Arg>
      static CODI_INLINE auto applyAdjointArg(Adjoint const& adjoint, Real const& result, Arg const& arg);

      // Interface implementation

      using Impl = CODI_DD(T_Impl, UnaryOperation);  ///< See UnaryOperation.

      /// \copydoc codi::ComputeOperation::applyTangent
      ///
      /// Forwards to applyTangentArg.
      template<size_t argNumber, typename Tangent, typename Arg>
      static CODI_INLINE auto applyTangent(Tangent const& tangent, Real const& result, Arg const& arg) {
        return Impl::applyTangentArg(tangent, result, arg);
      }

      /// \copydoc codi::ComputeOperation::applyAdjoint
      ///
      /// Forwards to applyAdjointArg.
      template<size_t argNumber, typename Adjoint, typename Arg>
      static CODI_INLINE auto applyAdjoint(Adjoint const& adjoint, Real const& result, Arg const& arg) {
        return Impl::applyAdjointArg(adjoint, result, arg);
      }
  };

  /**
   * @brief Implements UnaryOperation for functions where the gradient can be computed and transposed.
   *
   * Implementations need to define #primal and #gradient.
   *
   * See UnaryOperation for details.
   *
   * @tparam T_Real  Original primal value of the statement/expression.
   * @tparam T_Impl  Final implementation of the operation.
   */
  template<typename T_Real, typename T_Impl>
  struct UnaryJacobianOperation : UnaryOperation<T_Real, T_Impl> {
    public:

      using Real = CODI_DD(T_Real, double);  ///< See UnaryJacobianOperation.

      /// \copydoc codi::ComputeOperation::primal
      template<typename Arg>
      static CODI_INLINE Real primal(Arg const& arg);

      /// Compute \f$ \frac{\d f}{\d x} \f$.
      template<typename Arg>
      static CODI_INLINE auto gradient(Arg const& arg, Real const& result);

      // Interface implementation

      using Impl = CODI_DD(T_Impl, UnaryJacobianOperation);  ///< See UnaryOperation.

      /// \copydoc codi::UnaryOperation::applyTangentArg
      ///
      /// Calls gradient and multiplies with tangent.
      template<typename Tangent, typename Arg>
      static CODI_INLINE auto applyTangentArg(Tangent const& tangent, Real const& result, Arg const& arg) {
        return Impl::gradient(arg, result) * tangent;
      }

      /// \copydoc codi::UnaryOperation::applyAdjointArg
      ///
      /// Calls gradient, transposes the result and multiplies with adjoint.
      template<typename Adjoint, typename Arg>
      static CODI_INLINE Arg applyAdjointArg(Adjoint const& adjoint, Real const& result, Arg const& arg) {
        return ComputationTraits::transpose(Impl::gradient(arg, result)) * adjoint;
      }
  };

  /**
   * @brief Implements ComputeOperation for two arguments.
   *
   * Represents an operation \f$ w = f(a, b)\f$.
   *
   * Implementations need to define #primal, #applyTangentArgA, #applyTangentArgB, #applyAdjointArgA and
   * #applyAdjointArgB.
   *
   * See ComputeOperation for details.
   *
   * @tparam T_Real  Original primal value of the statement/expression.
   * @tparam T_Impl  Final implementation of the operation.
   */
  template<typename T_Real, typename T_Impl>
  struct BinaryOperation : public ComputeOperation<T_Real> {
    public:

      using Real = CODI_DD(T_Real, double);  ///< See BinaryOperation.

      /// \copydoc codi::ComputeOperation::primal
      template<typename ArgA, typename ArgB>
      static CODI_INLINE Real primal(ArgA const& argA, ArgB const& argB);

      /// Apply the forward AD mode with respect to argument a on the tangent and return the result.
      /// Computes \f$ \dot w = \frac{\d f}{\d a} \dot a \f$
      template<typename Tangent, typename ArgA, typename ArgB>
      static CODI_INLINE auto applyTangentArgA(Tangent const& tangent, Real const& result, ArgA const& argA,
                                               ArgB const& argB);

      /// Apply the reverse AD mode with respect to argument a on the adjoint and return the result.
      /// Computes \f$ \bar a = \frac{\d f}{\d a}^T \bar w \f$
      template<typename Adjoint, typename ArgA, typename ArgB>
      static CODI_INLINE auto applyAdjointArgA(Adjoint const& adjoint, Real const& result, ArgA const& argA,
                                               ArgB const& argB);

      /// Apply the forward AD mode with respect to argument b on the tangent and return the result.
      /// Computes \f$ \dot w = \frac{\d f}{\d b} \dot b \f$
      template<typename Tangent, typename ArgA, typename ArgB>
      static CODI_INLINE auto applyTangentArgB(Tangent const& tangent, Real const& result, ArgA const& argA,
                                               ArgB const& argB);

      /// Apply the reverse AD mode with respect to argument b on the adjoint and return the result.
      /// Computes \f$ \bar b = \frac{\d f}{\d b}^T \bar w \f$
      template<typename Adjoint, typename ArgA, typename ArgB>
      static CODI_INLINE auto applyAdjointArgB(Adjoint const& adjoint, Real const& result, ArgA const& argA,
                                               ArgB const& argB);

      // Interface implementation

      using Impl = CODI_DD(T_Impl, BinaryOperation);  ///< See BinaryOperation.

#ifndef DOXYGEN_DISABLE
      template<typename Func, bool derivB = false>
      struct CallSwitchTangent {
        public:
          template<typename... Args>
          CODI_INLINE static auto call(Args&&... args) {
            return Impl::applyTangentArgA(std::forward<Args>(args)...);
          }
      };

      template<typename Temp>
      struct CallSwitchTangent<Temp, true> {
        public:
          template<typename... Args>
          CODI_INLINE static auto call(Args&&... args) {
            return Impl::applyTangentArgB(std::forward<Args>(args)...);
          }
      };

      template<typename Func, bool derivB = false>
      struct CallSwitchAdjoint {
        public:
          template<typename... Args>
          CODI_INLINE static auto call(Args&&... args) {
            return Impl::applyAdjointArgA(std::forward<Args>(args)...);
          }
      };

      template<typename Temp>
      struct CallSwitchAdjoint<Temp, true> {
        public:
          template<typename... Args>
          CODI_INLINE static auto call(Args&&... args) {
            return Impl::applyAdjointArgB(std::forward<Args>(args)...);
          }
      };
#endif

      /// \copydoc codi::ComputeOperation::applyTangent
      ///
      /// Forwards to applyTangentArgA and applyTangentArgB.
      template<size_t argNumber, typename Tangent, typename ArgA, typename ArgB>
      static CODI_INLINE auto applyTangent(Tangent const& tangent, Real const& result, ArgA const& argA,
                                           ArgB const& argB) {
        return CallSwitchTangent<void, argNumber == 1>::call(tangent, result, argA, argB);
      }

      /// \copydoc codi::ComputeOperation::applyAdjoint
      ///
      /// Forwards to applyAdjointArgA and applyAdjointArgB.
      template<size_t argNumber, typename Adjoint, typename ArgA, typename ArgB>
      static CODI_INLINE auto applyAdjoint(Adjoint const& adjoint, Real const& result, ArgA const& argA,
                                           ArgB const& argB) {
        return CallSwitchAdjoint<void, argNumber == 1>::call(adjoint, result, argA, argB);
      }
  };

  /**
   * @brief Implements BinaryOperation for functions where the gradients can be computed and transposed.
   *
   * Implementations need to define #primal, #gradientA and #gradientB.
   *
   * See BinaryOperation for details.
   *
   * @tparam T_Real  Original primal value of the statement/expression.
   * @tparam T_Impl  Final implementation of the operation.
   */
  template<typename T_Real, typename T_Impl>
  struct BinaryJacobianOperation : public BinaryOperation<T_Real, T_Impl> {
    public:

      using Real = CODI_DD(T_Real, double);  ///< See BinaryJacobianOperation.

      /// \copydoc codi::ComputeOperation::primal
      template<typename ArgA, typename ArgB>
      static CODI_INLINE Real primal(ArgA const& argA, ArgB const& argB);

      /// Compute \f$ \frac{\d f}{\d a} \f$.
      template<typename ArgA, typename ArgB>
      static CODI_INLINE auto gradientA(ArgA const& argA, ArgB const& argB, Real const& result);

      /// Compute \f$ \frac{\d f}{\d b} \f$.
      template<typename ArgA, typename ArgB>
      static CODI_INLINE auto gradientB(ArgA const& argA, ArgB const& argB, Real const& result);

      // Interface implementation

      using Impl = CODI_DD(T_Impl, BinaryJacobianOperation);  ///< See BinaryJacobianOperation.

      /// \copydoc codi::BinaryOperation::applyTangentArgA
      ///
      /// Calls gradientA and multiplies with tangent.
      template<typename Tangent, typename ArgA, typename ArgB>
      static CODI_INLINE auto applyTangentArgA(Tangent const& tangent, Real const& result, ArgA const& argA,
                                               ArgB const& argB) {
        return Impl::gradientA(argA, argB, result) * tangent;
      }

      /// \copydoc codi::BinaryOperation::applyAdjointArgA
      ///
      /// Calls gradientB, transposes the result and multiplies with adjoint.
      template<typename Adjoint, typename ArgA, typename ArgB>
      static CODI_INLINE ArgA applyAdjointArgA(Adjoint const& adjoint, Real const& result, ArgA const& argA,
                                               ArgB const& argB) {
        return ComputationTraits::transpose(Impl::gradientA(argA, argB, result)) * adjoint;
      }

      /// \copydoc codi::BinaryOperation::applyTangentArgB
      ///
      /// Calls gradientB and multiplies with tangent.
      template<typename Tangent, typename ArgA, typename ArgB>
      static CODI_INLINE auto applyTangentArgB(Tangent const& tangent, Real const& result, ArgA const& argA,
                                               ArgB const& argB) {
        return Impl::gradientB(argA, argB, result) * tangent;
      }

      /// \copydoc codi::BinaryOperation::applyAdjointArgB
      ///
      /// Calls gradientB, transposes the result and multiplies with adjoint.
      template<typename Adjoint, typename ArgA, typename ArgB>
      static CODI_INLINE ArgB applyAdjointArgB(Adjoint const& adjoint, Real const& result, ArgA const& argA,
                                               ArgB const& argB) {
        return ComputationTraits::transpose(Impl::gradientB(argA, argB, result)) * adjoint;
      }
  };

  /**
   * @brief Represents an operator or function with an arbitrary number of arguments in the expression tree.
   *
   * See \ref Expressions "Expression" design documentation for details about the expression system in CoDiPack.
   *
   * @tparam T_Real  Original primal value of the statement/expression.
   * @tparam T_Operation  The logic for computing the primal value and Jacobians. Must implement BinaryOperation.
   * @tparam T_ArgExprs  The ExpressionInterface types of the arguments.
   */
  template<typename T_Real, template<typename> class T_Operation, typename... T_ArgExprs>
  struct ComputeExpression
      : public ExpressionInterface<T_Real, ComputeExpression<T_Real, T_Operation, T_ArgExprs...> > {
    public:
      using Real = CODI_DD(T_Real, double);  ///< See ComputeExpression.
      template<typename T>
      using Operation = CODI_DD(CODI_T(T_Operation<T>), CODI_T(ComputeOperation<T>));  ///< See ComputeExpression.
      using ArgExprs = std::tuple<T_ArgExprs...>;                                      ///< See ComputeExpression.

      using ArgReals = std::tuple<typename T_ArgExprs::Real...>;      ///< Real type of all arguments.
      using ArgStores = std::tuple<typename T_ArgExprs::StoreAs...>;  ///< Store type of all arguments.

      ArgStores args;  ///< Tuple of all expression arguments.

      Real result;  ///< Precomputed result.

      /// Constructor
      explicit ComputeExpression(ExpressionInterface<typename T_ArgExprs::Real, T_ArgExprs> const&... args)
          : args(args.cast()...), result(Operation<Real>::primal(args.cast().getValue()...)) {}

      /*******************************************************************************/
      /// @name Implementation of ExpressionInterface
      /// @{

      /// \copydoc codi::ComputeOperation::getMathRep
      CODI_INLINE std::string getMathRep() const {
        return Operation<Real>::getMathRep();
      }

      /// @}
      /*******************************************************************************/
      /// @name Implementation of ExpressionInterface
      /// @{

      using StoreAs = ComputeExpression;  ///< See ExpressionInterface.
      using ADLogic = typename ExpressionTraits::ValidateADLogic<
          typename T_ArgExprs::ADLogic...>::ADLogic;  ///< See ExpressionInterface.

      /// \copydoc codi::ExpressionInterface::getValue()
      CODI_INLINE Real const& getValue() const {
        return result;
      }

      /// \copydoc codi::ExpressionInterface::applyTangent()
      ///
      /// Forwards to the operator implementation.
      template<size_t argNumber, typename Tangent>
      CODI_INLINE Real applyTangent(Tangent const& tangent) const {
        return callTangent<argNumber>(tangent, result, args);
      }

      /// \copydoc codi::ExpressionInterface::applyAdjoint()
      ///
      /// Forwards to the operator implementation.
      template<size_t argNumber, typename Adjoint>
      CODI_INLINE auto applyAdjoint(Adjoint const& adjoint) const {
        return callAdjoint<argNumber>(adjoint, result, args);
      }

      /// @}
      /*******************************************************************************/
      /// @name Implementation of NodeInterface
      /// @{

      static size_t constexpr LinkCount = std::tuple_size<ArgReals>::value;  ///< See NodeInterface.

      /// \copydoc NodeInterface::getLink
      template<size_t argNumber>
      CODI_INLINE std::tuple_element_t<argNumber, ArgStores> const& getLink() const {
        return std::get<argNumber>(args);
      }
      /// @}

    private:

      // Unpack helpers

      template<size_t argNumber, typename EvalArg, typename Tuple, size_t... I>
      static CODI_INLINE auto callAdjoint(EvalArg const& evalArg, Real const& result, Tuple t,
                                          std::index_sequence<I...>) {
        return Operation<Real>::template applyAdjoint<argNumber>(evalArg, result, std::get<I>(t).getValue()...);
      }

      template<size_t argNumber, typename EvalArg, typename Tuple>
      static CODI_INLINE auto callAdjoint(EvalArg const& evalArg, Real const& result, Tuple const& t) {
        static constexpr auto size = std::tuple_size<Tuple>::value;
        return callAdjoint<argNumber>(evalArg, result, t, std::make_index_sequence<size>{});
      }

      template<size_t argNumber, typename EvalArg, typename Tuple, size_t... I>
      static CODI_INLINE auto callTangent(EvalArg const& evalArg, Real const& result, Tuple t,
                                          std::index_sequence<I...>) {
        return Operation<Real>::template applyTangent<argNumber>(evalArg, result, std::get<I>(t).getValue()...);
      }

      template<size_t argNumber, typename EvalArg, typename Tuple>
      static CODI_INLINE auto callTangent(EvalArg const& evalArg, Real const& result, Tuple const& t) {
        static constexpr auto size = std::tuple_size<Tuple>::value;
        return callTangent<argNumber>(evalArg, result, t, std::make_index_sequence<size>{});
      }
  };
}
