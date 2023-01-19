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

#include <complex>

#include "../../config.h"
#include "../../misc/macros.hpp"
#include "../../traits/expressionTraits.hpp"
#include "../unaryExpression.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  struct ReduceToReal {};  ///< Placeholder to identify the operation on the Jacobian.

  /**
   * @brief Returns a proxy object for the gradient that implements the operation in the multiplication of the proxy.
   *
   * In CoDiPack the Jacobian is applied to other Jacobians via multiplication. Since the adjoint of this operation can
   * not be described by a complex or floating point value, that would apply the same logic, a place holder is required.
   * With this placeholder the multiplication operation is specialized and the custom logic is evaluated.
   *
   * @tparam T_Real  Original primal value of the statement/expression. (E.g. for double for complex<double>.)
   */
  template<typename T_Real>
  struct OperationAdjointComplexToRealCast : public UnaryOperation<T_Real> {
    public:

      using Real = CODI_DD(T_Real, double);  ///< See OperationAdjointComplexToRealCast.
      using Jacobian = ReduceToReal;         ///< See OperationAdjointComplexToRealCast.

      /// \copydoc codi::UnaryOperation::primal
      template<typename Arg>
      static CODI_INLINE Real primal(Arg const& arg) {
        return arg;
      }

      /// See OperationAdjointComplexToRealCast.
      template<typename Arg>
      static CODI_INLINE ReduceToReal gradient(Arg const& arg, Real const& result) {
        CODI_UNUSED(arg, result);

        return ReduceToReal{};
      }
  };

  /// Expression that converts in the adjoint evaluation a complex to the real part. See
  /// OperationAdjointComplexToRealCast for details.
  template<typename T_Real, typename T_Arg>
  using AdjointComplexToRealCast = UnaryExpression<T_Real, T_Arg, OperationAdjointComplexToRealCast>;

  /// See codi::OperationAdjointComplexToRealCast.
  template<typename Real>
  CODI_INLINE Real operator*(ReduceToReal, std::complex<Real> const& adjoint) {
    return adjoint.real();
  }

  /// See codi::OperationAdjointComplexToRealCast.
  template<typename Real, typename Arg>
  CODI_INLINE ExpressionTraits::ActiveResult<Real, typename Arg::ADLogic> operator*(
      ReduceToReal, ExpressionInterface<std::complex<Real>, Arg> const& adjoint) {
    return adjoint.cast().real();
  }
}
