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
 * In order to include this file the user has to define the preprocessor macro OPERATOR.
 * OPERATOR contains the name of the comparison operator without the 'operator' classifier.
 * e.g. '<=' or '>'.
 *
 * The define OPERATOR will be undefined at the end of this template.
 */

#ifndef OPERATOR
  #error Please define the name of the operator.
#endif

// Create a correct include environment for viewing and programming in an IDE.
#ifndef OPERATOR
  #define PROXY_OUTER

  #include "../../config.h"
  #include "../../misc/macros.hpp"
  #include "../expressionInterface.hpp"
  #define OPERATOR ==

namespace codi {
#endif

  /// Function overload for operator OPERATOR.
  template<typename Real, typename ArgA, typename ArgB>
  CODI_INLINE bool operator OPERATOR(ExpressionInterface<Real, ArgA> const& argA,
                                     ExpressionInterface<Real, ArgB> const& argB) {
    return RealTraits::getPassiveValue(argA.cast()) OPERATOR RealTraits::getPassiveValue(argB.cast());
  }

#define PASSIVE_TYPE RealTraits::PassiveReal<Real>
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

// Create a correct include environment for viewing and programming in an IDE.
#ifdef PROXY_OUTER
  #undef PROXY_OUTER
}
#endif

#undef OPERATOR
