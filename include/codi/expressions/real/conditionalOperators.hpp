/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2025 Chair for Scientific Computing (SciComp), University of Kaiserslautern-Landau
 * Homepage: http://scicomp.rptu.de
 * Contact:  Prof. Nicolas R. Gauger (codi@scicomp.uni-kl.de)
 *
 * Lead developers: Max Sagebaum, Johannes Blühdorn (SciComp, University of Kaiserslautern-Landau)
 *
 * This file is part of CoDiPack (http://scicomp.rptu.de/software/codi).
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
#pragma once

#include "../../config.h"
#include "../../misc/exceptions.hpp"
#include "../../misc/macros.hpp"
#include "../../traits/realTraits.hpp"
#include "../expressionInterface.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /*******************************************************************************/
  /// @name Builtin binary comparison operators
  /// @{

#define RETURN bool
#define OPERATOR ==
#include "conditionalBinaryOverloads.tpp"

#define RETURN bool
#define OPERATOR !=
#include "conditionalBinaryOverloads.tpp"

#define RETURN bool
#define OPERATOR >
#include "conditionalBinaryOverloads.tpp"

#define RETURN bool
#define OPERATOR <
#include "conditionalBinaryOverloads.tpp"

#define RETURN bool
#define OPERATOR >=
#include "conditionalBinaryOverloads.tpp"

#define RETURN bool
#define OPERATOR <=
#include "conditionalBinaryOverloads.tpp"

#define RETURN bool
#define OPERATOR &&
#include "conditionalBinaryOverloads.tpp"

#define RETURN bool
#define OPERATOR ||
#include "conditionalBinaryOverloads.tpp"

#if CODI_HasCpp20
  #define RETURN std::partial_ordering
  #define OPERATOR <=>
  #include "conditionalBinaryOverloads.tpp"
#endif

  /// @}
  /*******************************************************************************/
  /// @name Builtin unary comparison operators
  /// @{

#define OPERATOR !
#include "conditionalUnaryOverloads.tpp"

  /// @}
}
