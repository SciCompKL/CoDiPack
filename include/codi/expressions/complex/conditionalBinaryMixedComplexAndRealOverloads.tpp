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

#ifndef OPERATOR
  #error Please define the operator for the comparison.
#endif

// Create a correct include environment for viewing and programming in an IDE
#ifndef OPERATOR
  #include <complex>

  #include "../../config.h"
  #include "../../misc/macros.hpp"
  #include "../../traits/realTraits.hpp"
  #include "../expressionInterface.hpp"

  #define OPERATOR ==

namespace codi {
#endif

  // Do not need to define complex complex binding, they are handled by the default real definitions

  /// Function overload for OPERATOR(complex, const real).
  template<typename Real, typename ArgA>
  CODI_INLINE bool operator OPERATOR(ExpressionInterface<std::complex<Real>, ArgA> const& argA,
                                     RealTraits::PassiveReal<Real> const& argB) {
    return RealTraits::getPassiveValue(argA.cast()) OPERATOR argB;
  }

  /// Function overload for OPERATOR(const real, complex).
  template<typename Real, typename ArgB>
  CODI_INLINE bool operator OPERATOR(RealTraits::PassiveReal<Real> const& argA,
                                     ExpressionInterface<std::complex<Real>, ArgB> const& argB) {
    return argA OPERATOR RealTraits::getPassiveValue(argB.cast());
  }

// Create a correct include environment for viewing and programming in an IDE
#ifndef OPERATOR
}
#endif

#undef OPERATOR
