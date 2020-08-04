#pragma once

#include "../../aux/exceptions.hpp"
#include "../../aux/macros.hpp"
#include "../../config.h"
#include "../../traits/realTraits.hpp"
#include "../expressionInterface.hpp"

/** \copydoc codi::Namespace */
namespace codi {


  /*******************************************************************************
   * Section: Standard binary comparison operators
   *
   * Description: TODO
   *
   */

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

  /*******************************************************************************
   * Section: Standard prefix unary comparison operators
   *
   * Description: TODO
   *
   */

  #define OPERATOR !
  #include "conditionalUnaryOverloads.tpp"
}

