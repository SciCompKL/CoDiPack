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

/** \copydoc codi::Namespace */
namespace codi {

  /** Interface for injecting member operations in expressions
   *
   *  The ExpressionInterface extends from this interface and therefore all methods defined in specializations of
   *  this class are available in expressions with the specialized Real type.
   *
   *  E.g. the complex type has the member operations real() and imag(). Which makes the following code possible.
   *  \code
   *  std::complex<double> a,b;
   *  double z = (a + b).real();
   *  \endcode
   *  In an expression framework, the member operator could be implemented in the AggregatedActiveType specialization
   *  for std::complex but this will not make the member operator available in the expression (a+b). The compiler would
   *  throw an error that the member operator `real()` is not available in a BinaryExpression.
   *
   *  By the specialization of this class, member methods can be injected into arbitrary expression implementations.
   *  This includes e.g. ActiveType, AggregatedActiveType, BinaryExpression, etc. In the complex case the specialization
   *  can be viewed in stdComplex.hpp.
   *
   * @tparam T_Real  Original primal value of the statement/expression.
   * @tparam T_Impl  Class implementing this interface.
   */
  template<typename T_Real, typename T_Impl, typename = void>
  struct ExpressionMemberOperations {
      using Real = CODI_DD(T_Real, double);    ///< See ExpressionMemberOperations.
      using Impl = CODI_DD(T_Impl, CODI_ANY);  ///< See ExpressionMemberOperations.

    protected:

      /// Cast to implementation.
      CODI_INLINE Impl const& cast() const {
        return static_cast<Impl const&>(*this);
      }
  };
}
