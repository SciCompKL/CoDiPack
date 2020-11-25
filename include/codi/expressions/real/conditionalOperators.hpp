#pragma once

#include "../../aux/exceptions.hpp"
#include "../../aux/macros.hpp"
#include "../../config.h"
#include "../../traits/realTraits.hpp"
#include "../expressionInterface.hpp"

/** \copydoc codi::Namespace */
namespace codi {


  /*******************************************************************************/
  /// @name Builtin binary comparison operators
  /// @{

  #define OPERATOR ==
  #include "conditionalBinaryOverloads.tpp"

  #define OPERATOR !=
  #include "conditionalBinaryOverloads.tpp"

  #define OPERATOR >
  #include "conditionalBinaryOverloads.tpp"

  #define OPERATOR <
  #include "conditionalBinaryOverloads.tpp"

  #define OPERATOR >=
  #include "conditionalBinaryOverloads.tpp"

  #define OPERATOR <=
  #include "conditionalBinaryOverloads.tpp"

  #define OPERATOR &&
  #include "conditionalBinaryOverloads.tpp"

  #define OPERATOR ||
  #include "conditionalBinaryOverloads.tpp"

  /// @}
  /*******************************************************************************/
  /// @name Builtin unary comparison operators
  /// @{

  #define OPERATOR !
  #include "conditionalUnaryOverloads.tpp"

  /// @}
}

