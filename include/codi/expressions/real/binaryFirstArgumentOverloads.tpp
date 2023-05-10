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
#ifndef SECOND_ARG_TYPE
  #error Please define the type of the second argument.
#endif
#ifndef SECOND_ARG_CONVERSION
  #error Please define the conversion operations for the second argument.
#endif

// Create a correct include environment for viewing and programming in an IDE.
#ifndef FUNCTION
  #define PROXY

  #include "../../config.h"
  #include "../../misc/macros.hpp"
  #include "../../traits/realTraits.hpp"
  #include "../binaryExpression.hpp"
  #include "../constantExpression.hpp"
  #include "../expressionInterface.hpp"
  #define OPERATION_LOGIC BinaryOperation
  #define FUNCTION func
  #define SECOND_ARG_TYPE double
  #define SECOND_ARG_CONVERSION ConstantDataConversion

namespace codi {
#endif

  /// Function overload for FUNCTION.
  template<typename Real, typename ArgA>
  CODI_INLINE BinaryExpression<Real, ArgA, ConstantExpression<SECOND_ARG_TYPE, SECOND_ARG_CONVERSION>, OPERATION_LOGIC>
  FUNCTION(ExpressionInterface<Real, ArgA> const& argA, SECOND_ARG_TYPE const& argB) {
    return BinaryExpression<Real, ArgA, ConstantExpression<SECOND_ARG_TYPE, SECOND_ARG_CONVERSION>, OPERATION_LOGIC>(
        argA, ConstantExpression<SECOND_ARG_TYPE, SECOND_ARG_CONVERSION>(argB));
  }

// Create a correct include environment for viewing and programming in an IDE.
#ifdef PROXY
  #undef PROXY
}
#endif

#undef FUNCTION
#undef OPERATION_LOGIC
#undef SECOND_ARG_TYPE
#undef SECOND_ARG_CONVERSION
