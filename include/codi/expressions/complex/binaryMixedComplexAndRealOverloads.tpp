/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2026 Chair for Scientific Computing (SciComp), University of Kaiserslautern-Landau
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

// clang-format off
/*
 * This file defines function overloads of the kind:
 *
 * The (complex, complex) case is not defined, since it is handled by the real definitions.
 * ComputeExpression<std::complex<Real>, ...> FUNCTION(ExpressionInterface<std::complex<Real>, ArgA> const&, ExpressionInterface<Real, ArgB> const&);
 * ComputeExpression<std::complex<Real>, ...> FUNCTION(ExpressionInterface<std::complex<Real>, ArgA> const&, PassiveReal const&);
 * ComputeExpression<std::complex<Real>, ...> FUNCTION(std::complex<PassiveReal> const&,                     ExpressionInterface<Real, ArgB> const&);
 * ComputeExpression<std::complex<Real>, ...> FUNCTION(ExpressionInterface<Real, ArgA> const&,               ExpressionInterface<std::complex<Real>, ArgB> const&);
 * ComputeExpression<std::complex<Real>, ...> FUNCTION(ExpressionInterface<Real, ArgA> const&,               std::complex<PassiveReal> const&);
 * ComputeExpression<std::complex<Real>, ...> FUNCTION(PassiveReal const&,                                   ExpressionInterface<std::complex<Real>, ArgB> const&);
 *
 * In order to include this file the user has to define the preprocessor macros OPERATION_LOGIC and FUNCTION.
 * OPERATION_LOGIC contains the name of the operation logic class. FUNCTION represents the normal name of that function
 * e.g. 'operator -' or 'sin'.
 *
 * The defines OPERATION_LOGIC and FUNCTION will be undefined at the end of this template.
 *
 * Prior to including this file, the user has to implement the operation's primal and derivative logic according to
 * BinaryOpInterface.
 */
// clang-format on

#ifndef OPERATION_LOGIC
  #error Please define a name for the binary expression.
#endif
#ifndef FUNCTION
  #error Please define the primal function representation.
#endif

// Create a correct include environment for viewing and programming in an IDE.
#ifndef FUNCTION
  #include <complex>

  #include "../../config.h"
  #include "../../misc/macros.hpp"
  #include "../../traits/realTraits.hpp"
  #include "../activeType.hpp"
  #include "../constantExpression.hpp"
  #include "../expressionInterface.hpp"
  #include "realToComplexCast.hpp"

  #define OPERATION_LOGIC BinaryJacobianOperation
  #define FUNCTION func

namespace codi {
#endif

  // No need to define (complex, complex) bindings, they are handled by the default real definitions.

  // Define (complex, real) bindings

  /// Function overload for FUNCTION(complex, real).
  template<typename Real, typename ArgA, typename ArgB>
  CODI_INLINE auto FUNCTION(ExpressionInterface<std::complex<Real>, ArgA> const& argA,
                            ExpressionInterface<Real, ArgB> const& argB) {
    return ComputeExpression<std::complex<Real>, OPERATION_LOGIC, ArgA, RealToComplexCast<Real, ArgB>>(
        argA, RealToComplexCast<Real, ArgB>(argB));
  }

  /// Function overload for FUNCTION(complex, passive real).
  template<typename Real, typename ArgA>
  CODI_INLINE auto FUNCTION(ExpressionInterface<std::complex<Real>, ArgA> const& argA,
                            RealTraits::PassiveReal<Real> const& argB) {
    return ComputeExpression<std::complex<Real>, OPERATION_LOGIC, ArgA,
                             ConstantExpression<RealTraits::PassiveReal<Real>>>(
        argA, ConstantExpression<RealTraits::PassiveReal<Real>>(argB));
  }

  /// Function overload for FUNCTION(passive complex, real).
  template<typename Real, typename ArgB>
  CODI_INLINE auto FUNCTION(std::complex<RealTraits::PassiveReal<Real>> const& argA,
                            ExpressionInterface<Real, ArgB> const& argB) {
    return ComputeExpression<std::complex<Real>, OPERATION_LOGIC,
                             ConstantExpression<std::complex<RealTraits::PassiveReal<Real>>>,
                             RealToComplexCast<Real, ArgB>>(
        ConstantExpression<std::complex<RealTraits::PassiveReal<Real>>>(argA), RealToComplexCast<Real, ArgB>(argB));
  }

  // Define (real, complex) bindings

  /// Function overload for FUNCTION(real, complex).
  template<typename Real, typename ArgA, typename ArgB>
  CODI_INLINE auto FUNCTION(ExpressionInterface<Real, ArgA> const& argA,
                            ExpressionInterface<std::complex<Real>, ArgB> const& argB) {
    return ComputeExpression<std::complex<Real>, OPERATION_LOGIC, RealToComplexCast<Real, ArgA>, ArgB>(
        RealToComplexCast<Real, ArgA>(argA), argB);
  }

  /// Function overload for FUNCTION(real, passive complex).
  template<typename Real, typename ArgA>
  CODI_INLINE auto FUNCTION(ExpressionInterface<Real, ArgA> const& argA,
                            std::complex<RealTraits::PassiveReal<Real>> const& argB) {
    return ComputeExpression<std::complex<Real>, OPERATION_LOGIC, RealToComplexCast<Real, ArgA>,
                             ConstantExpression<std::complex<RealTraits::PassiveReal<Real>>>>(
        RealToComplexCast<Real, ArgA>(argA), ConstantExpression<std::complex<RealTraits::PassiveReal<Real>>>(argB));
  }

  /// Function overload for FUNCTION(passive real, complex).
  template<typename Real, typename ArgB>
  CODI_INLINE auto FUNCTION(RealTraits::PassiveReal<Real> const& argA,
                            ExpressionInterface<std::complex<Real>, ArgB> const& argB) {
    return ComputeExpression<std::complex<Real>, OPERATION_LOGIC, ConstantExpression<RealTraits::PassiveReal<Real>>,
                             ArgB>(ConstantExpression<RealTraits::PassiveReal<Real>>(argA), argB);
  }

// Create a correct include environment for viewing and programming in an IDE.
#ifndef FUNCTION
}
#endif

#undef FUNCTION
#undef OPERATION_LOGIC
