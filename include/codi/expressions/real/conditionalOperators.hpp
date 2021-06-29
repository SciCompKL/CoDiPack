/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2021 Chair for Scientific Computing (SciComp), TU Kaiserslautern
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
