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

#include <iostream>

#include "../config.h"
#include "../misc/macros.hpp"
#include "../traits/expressionTraits.hpp"
#include "../traits/realTraits.hpp"
#include "expressionMemberOperations.hpp"
#include "logic/nodeInterface.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Base class for all CoDiPack expressions.
   *
   * See \ref Expressions "Expression" design documentation for details about the expression system in CoDiPack.
   *
   * This interface resembles a rvalue in C++.
   *
   * @tparam T_Real  Original primal value of the statement/expression.
   * @tparam T_Impl  Class implementing this interface.
   */
  template<typename T_Real, typename T_Impl>
  struct ExpressionInterface : public NodeInterface<T_Impl>, public ExpressionMemberOperations<T_Real, T_Impl> {
    public:

      using Real = CODI_DD(T_Real, double);               ///< See ExpressionInterface.
      using Impl = CODI_DD(T_Impl, ExpressionInterface);  ///< See ExpressionInterface.

      /// AD logic that governs the expression. Needs to be the same for all inputs of the expression.
      using ADLogic = CODI_UNDEFINED;

      /// Constructor
      ExpressionInterface() = default;

      /// Constructor
      ExpressionInterface(ExpressionInterface const&) = default;

      /// Cast to the implementation.
      CODI_INLINE Impl const& cast() const {
        return static_cast<Impl const&>(*this);
      }

#if CODI_ImplicitConversion
      /// Implicit cast for CoDiPack expressions.
      CODI_INLINE operator const Real() const {
        Warning::implicitCast<Config::ImplicitConversionWarning>();

        return cast().getValue();
      }
#endif

      /*******************************************************************************/
      /// @name Interface definition
      /// @{

      using StoreAs = ExpressionInterface;  ///< Defines how this expression is stored in an expression tree.

      /// Compute the primal value that is usually evaluated by the statement/expression.
      CODI_INLINE Real const getValue() const;

      /** Apply the AD forward mode on the expression with respect to the given parameter.
       *
       *  This is just the local forward mode application and not the one for the whole expression tree.
       *
       *  Does not need to be implemented for expressions with \c NodeInterface::LinkCount = 0 \c .
       *
       *  @return The type of the result or a compatible vector type. E.g. Real or Direction<Real>.
       *
       *  @tparam Tangent  The type is the Real type of the selected argument or a compatible vector type. E.g. for
       *  \c Real f(complex<Real>, Real) \c the type with \c argNumber=0 \c is \c complex<Real> \c or
       *  \c Direction<complex<Real>> \c , with \c argNumber=1 \c it is \c Real \c or \c Direction<Real> \c .
       */
      template<size_t argNumber, typename Tangent>
      CODI_INLINE auto applyTangent(Tangent const& tangent) const;

      /**  Apply the AD reverse mode on the expression with respect to the given parameter.
       *
       *  This is just the local reverse mode application and not the one for the whole expression tree.
       *
       *  Does not need to be implemented for expressions with \c NodeInterface::LinkCount = 0 \c .
       *
       *  @return The type is the Real type of the selected argument or a compatible vector type. E.g. for
       *  \c Real f(complex<Real>, Real) \c the type with \c argNumber=0 \c is \c complex<Real> \c or
       *  \c Direction<complex<Real>> \c , with \c argNumber=1 \c it is \c Real \c or \c Direction<Real> \c .
       *
       *  @tparam Adjoint  The type of the result or a compatible vector type. E.g. Real or Direction<Real>.
       */
      template<size_t argNumber, typename Adjoint>
      CODI_INLINE auto applyAdjoint(Adjoint const& adjoint) const;

      /// @}

    private:
      ExpressionInterface& operator=(ExpressionInterface const&) = delete;
  };

#ifndef DOXYGEN_DISABLE
  template<typename T_Type>
  struct RealTraits::TraitsImplementation<T_Type, ExpressionTraits::EnableIfExpression<T_Type>> {
    public:

      using Type = CODI_DD(T_Type, CODI_T(ExpressionInterface<double, T_Type>));
      using Real = typename Type::Real;

      using PassiveReal = RealTraits::PassiveReal<Real>;

      static int constexpr MaxDerivativeOrder = 1 + RealTraits::MaxDerivativeOrder<Real>();

      static CODI_INLINE PassiveReal getPassiveValue(Type const& v) {
        return RealTraits::getPassiveValue(v.getValue());
      }
  };
#endif

  /// Write the primal value to the stream.
  template<typename Expr>
  ExpressionTraits::EnableIfExpression<Expr, std::ostream>& operator<<(std::ostream& out, Expr const& v) {
    out << v.getValue();

    return out;
  }
}
