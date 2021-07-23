/*
 * In order to include this file the user has to define the preprocessor macro OPERATION_LOGIC and FUNCTION.
 * OPERATION_LOGIC contains the name of the operation logic class. FUNCTION represents the normal name of that function
 * e.g. 'operator -' or 'sin'.
 *
 * The defines OPERATION_LOGIC and FUNCTION will be undefined at the end of this template.
 *
 * Prior to including this file, the user has to implement the operation's primal and derivative logic according to
 * UnaryOpInterface.
 */

#ifndef OPERATION_LOGIC
  #error Please define a name for the unary expression.
#endif
#ifndef FUNCTION
  #error Please define the primal function representation.
#endif

// Create a correct include environment for viewing and programming in an IDE.
#ifndef FUNCTION
  #define PROXY

  #include "../../aux/macros.hpp"
  #include "../../config.h"
  #include "../expressionInterface.hpp"
  #include "../unaryExpression.hpp"
  #define FUNCTION func
  #define OPERATION_LOGIC UnaryOperation

namespace codi {
#endif

  /// Function overload for FUNCTION.
  template<typename Real, typename Arg>
  CODI_INLINE UnaryExpression<Real, Arg, OPERATION_LOGIC> FUNCTION(ExpressionInterface<Real, Arg> const& arg) {
    return UnaryExpression<Real, Arg, OPERATION_LOGIC>(arg);
  }

// Create a correct include environment for viewing and programming in an IDE.
#ifdef PROXY
  #undef PROXY
}
#endif

#undef FUNCTION
#undef OPERATION_LOGIC
