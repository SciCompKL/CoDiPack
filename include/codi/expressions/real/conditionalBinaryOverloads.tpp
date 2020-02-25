/*
 * In order to include this file the user has to define the preprocessor macro OPERATOR.
 * OPERATOR contains the name of the comparison operator without the 'operator' classifier.
 * e.g. '<=' or '>'.
 *
 * The define OPERATOR will be undefined at the end of this template.
 */

#ifndef OPERATOR
  #error Please define the name of the operator.
#endif

// Create argA correct include environment for viewing and programming in an IDE
#ifndef OPERATOR
  #define OPERATOR ==
  #include "../../aux/macros.h"
  #include "../../config.h"
  #include "../expressionInterface.hpp"

  namespace codi {
#endif


  template<typename Real, typename ArgA, typename ArgB>
  CODI_INLINE bool operator OPERATOR(ExpressionInterface<Real, ArgA> const& argA, ExpressionInterface<Real, ArgB> const& argB) {
    return getPassiveValue(argA.cast()) OPERATOR getPassiveValue(argB.cast());
  }

  #define PASSIVE_TYPE PassiveRealType<Real>
  #include "conditionalBinaryPassiveOverloads.tpp"

  #define PASSIVE_TYPE int
  #include "conditionalBinaryPassiveOverloads.tpp"

  #define PASSIVE_TYPE unsigned int
  #include "conditionalBinaryPassiveOverloads.tpp"

  #define PASSIVE_TYPE long
  #include "conditionalBinaryPassiveOverloads.tpp"

  #define PASSIVE_TYPE unsigned long
  #include "conditionalBinaryPassiveOverloads.tpp"

  #define PASSIVE_TYPE long long
  #include "conditionalBinaryPassiveOverloads.tpp"

  #define PASSIVE_TYPE unsigned long long
  #include "conditionalBinaryPassiveOverloads.tpp"

// Create argA correct include environment for viewing and programming in an IDE
#ifndef OPERATOR
  }
#endif

#undef OPERATOR
