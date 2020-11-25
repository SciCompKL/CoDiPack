/*
 * This file should only be used in conditionalBinaryOverloads.tpp. It uses the defintions of this file.
 * In addition PASSIVE_TYPE needs to be declared. It defines the passive type for which the operator is overloaded.
 *
 * The define PASSIVE_TYPE will be undefined at the end of this template.
 */

#ifndef OPERATOR
  #error Please define the name of the operator.
#endif

#ifndef PASSIVE_TYPE
  #error Please define the passive type for the overloads.
#endif

// Create a correct include environment for viewing and programming in an IDE
#ifndef OPERATOR
  #define PROXY

  #include "../../aux/macros.hpp"
  #include "../../config.h"
  #include "../../traits/realTraits.hpp"
  #include "../expressionInterface.hpp"
  #define OPERATOR ==
  #define PASSIVE_TYPE double

  namespace codi {
#endif

    /// Function overload for operator OPERATOR
  template<typename Real, typename ArgA>
  CODI_INLINE bool operator OPERATOR(ExpressionInterface<Real, ArgA> const& argA, PASSIVE_TYPE const& argB) {
    return RealTraits::getPassiveValue(argA.cast()) OPERATOR argB;
  }

  /// Function overload for operator OPERATOR
  template<typename Real, typename ArgB>
  CODI_INLINE bool operator OPERATOR(PASSIVE_TYPE const& argA, ExpressionInterface<Real, ArgB> const& argB) {
    return argA OPERATOR RealTraits::getPassiveValue(argB.cast());
  }

// Create a correct include environment for viewing and programming in an IDE
#ifdef PROXY
  #undef PROXY
  }
#endif

#undef PASSIVE_TYPE
