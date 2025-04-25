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
#include "../computeExpression.hpp"

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
  template<typename T_ComplexReal>
  struct OperationRealToComplexCast : public UnaryOperation<T_ComplexReal, OperationRealToComplexCast<T_ComplexReal>> {
    public:

      using ComplexReal = CODI_DD(T_ComplexReal, std::complex<double>);  ///< See OperationRealToComplexCast.
      using Jacobian = ReduceToReal;                                     ///< See OperationRealToComplexCast.

      using Real = typename ComplexReal::value_type;  ///< Inner type of the complex value.

      /// \copydoc codi::UnaryOperation::primal
      template<typename Arg>
      static CODI_INLINE ComplexReal primal(Arg const& arg) {
        return arg;
      }

      /// \copydoc UnaryOperation::applyTangentArg
      template<typename Tangent, typename Arg>
      static CODI_INLINE ComplexReal applyTangentArg(Tangent const& tangent, ComplexReal const& result,
                                                     Arg const& arg) {
        CODI_UNUSED(result, arg);
        return tangent;
      }

      /// \copydoc UnaryOperation::applyAdjointArg
      template<typename Adjoint, typename Arg>
      static CODI_INLINE Real applyAdjointArg(Adjoint const& adjoint, ComplexReal const& result, Arg const& arg) {
        CODI_UNUSED(result, arg);

        return adjoint.real();
      }


      /// \copydoc codi::BinaryJacobianOperation::getMathRep()
      static CODI_INLINE std::string getMathRep() {
        return "()";
      }
  };

  /// Expression that converts in the adjoint evaluation a complex to the real part. See
  /// OperationRealToComplexCast for details.
  template<typename T_Real, typename T_Arg>
  using RealToComplexCast = ComputeExpression<std::complex<T_Real>, OperationRealToComplexCast, T_Arg>;
}
