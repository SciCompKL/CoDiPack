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

/*
 * In order to include this file the user has to define the preprocessor macros OPERATION_LOGIC and FUNCTION.
 * OPERATION_LOGIC contains the name of the operation logic class. FUNCTION represents the normal name of that function
 * e.g. 'operator -' or 'sin'.
 *
 * The defines OPERATION_LOGIC and FUNCTION will be undefined at the end of this template.
 *
 * Prior to including this file, the user has to implement the operation's primal and derivative logic according to
 * BinaryOpInterface.
 */

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
  #include "../binaryExpression.hpp"
  #include "../constantExpression.hpp"
  #include "../expressionInterface.hpp"
  #include "adjointComplexToRealCast.hpp"

  #define OPERATION_LOGIC BinaryOperation
  #define FUNCTION func

namespace codi {
#endif

  // No need to define complex complex bindings, they are handled by the default real definitions.

  // Define complex real bindings

  /// Function overload for FUNCTION(complex, real).
  template<typename Real, typename ArgA, typename ArgB>
  CODI_INLINE BinaryExpression<std::complex<Real>, ArgA, AdjointComplexToRealCast<Real, ArgB>, OPERATION_LOGIC>
  FUNCTION(ExpressionInterface<std::complex<Real>, ArgA> const& argA, ExpressionInterface<Real, ArgB> const& argB) {
    return BinaryExpression<std::complex<Real>, ArgA, AdjointComplexToRealCast<Real, ArgB>, OPERATION_LOGIC>(
        argA, AdjointComplexToRealCast<Real, ArgB>(argB));
  }

  /// Function overload for FUNCTION(const complex, real).
  template<typename Real, typename ArgA>
  CODI_INLINE
      BinaryExpression<std::complex<Real>, ArgA, ConstantExpression<RealTraits::PassiveReal<Real>>, OPERATION_LOGIC>
      FUNCTION(ExpressionInterface<std::complex<Real>, ArgA> const& argA, RealTraits::PassiveReal<Real> const& argB) {
    return BinaryExpression<std::complex<Real>, ArgA, ConstantExpression<RealTraits::PassiveReal<Real>>,
                            OPERATION_LOGIC>(argA, ConstantExpression<RealTraits::PassiveReal<Real>>(argB));
  }

  /// Function overload for FUNCTION(complex, const real).
  template<typename Real, typename ArgB>
  CODI_INLINE BinaryExpression<std::complex<Real>, ConstantExpression<std::complex<RealTraits::PassiveReal<Real>>>,
                               AdjointComplexToRealCast<Real, ArgB>, OPERATION_LOGIC>
  FUNCTION(std::complex<RealTraits::PassiveReal<Real>> const& argA, ExpressionInterface<Real, ArgB> const& argB) {
    return BinaryExpression<std::complex<Real>, ConstantExpression<std::complex<RealTraits::PassiveReal<Real>>>,
                            AdjointComplexToRealCast<Real, ArgB>, OPERATION_LOGIC>(
        ConstantExpression<std::complex<RealTraits::PassiveReal<Real>>>(argA),
        AdjointComplexToRealCast<Real, ArgB>(argB));
  }

  // Define real complex bindings

  /// Function overload for FUNCTION(real, complex).
  template<typename Real, typename ArgA, typename ArgB>
  CODI_INLINE BinaryExpression<std::complex<Real>, AdjointComplexToRealCast<Real, ArgA>, ArgB, OPERATION_LOGIC>
  FUNCTION(ExpressionInterface<Real, ArgA> const& argA, ExpressionInterface<std::complex<Real>, ArgB> const& argB) {
    return BinaryExpression<std::complex<Real>, AdjointComplexToRealCast<Real, ArgA>, ArgB, OPERATION_LOGIC>(
        AdjointComplexToRealCast<Real, ArgA>(argA), argB);
  }

  /// Function overload for FUNCTION(real, const complex).
  template<typename Real, typename ArgA>
  CODI_INLINE BinaryExpression<std::complex<Real>, AdjointComplexToRealCast<Real, ArgA>,
                               ConstantExpression<std::complex<RealTraits::PassiveReal<Real>>>, OPERATION_LOGIC>
  FUNCTION(ExpressionInterface<Real, ArgA> const& argA, std::complex<RealTraits::PassiveReal<Real>> const& argB) {
    return BinaryExpression<std::complex<Real>, AdjointComplexToRealCast<Real, ArgA>,
                            ConstantExpression<std::complex<RealTraits::PassiveReal<Real>>>, OPERATION_LOGIC>(
        AdjointComplexToRealCast<Real, ArgA>(argA),
        ConstantExpression<std::complex<RealTraits::PassiveReal<Real>>>(argB));
  }

  /// Function overload for FUNCTION(const real, complex).
  template<typename Real, typename ArgB>
  CODI_INLINE
      BinaryExpression<std::complex<Real>, ConstantExpression<RealTraits::PassiveReal<Real>>, ArgB, OPERATION_LOGIC>
      FUNCTION(RealTraits::PassiveReal<Real> const& argA, ExpressionInterface<std::complex<Real>, ArgB> const& argB) {
    return BinaryExpression<std::complex<Real>, ConstantExpression<RealTraits::PassiveReal<Real>>, ArgB,
                            OPERATION_LOGIC>(ConstantExpression<RealTraits::PassiveReal<Real>>(argA), argB);
  }

// Create a correct include environment for viewing and programming in an IDE.
#ifndef FUNCTION
}
#endif

#undef FUNCTION
#undef OPERATION_LOGIC
