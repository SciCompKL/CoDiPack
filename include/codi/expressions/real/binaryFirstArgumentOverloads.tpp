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

  #include "../../aux/macros.hpp"
  #include "../../config.h"
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
  CODI_INLINE BinaryExpression<Real, ArgA, ConstantExpression<SECOND_ARG_TYPE, SECOND_ARG_CONVERSION>, OPERATION_LOGIC> FUNCTION(
      ExpressionInterface<Real, ArgA> const& argA, SECOND_ARG_TYPE const& argB) {
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
