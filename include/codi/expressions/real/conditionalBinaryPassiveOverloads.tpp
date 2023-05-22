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

// Create a correct include environment for viewing and programming in an IDE.
#ifndef OPERATOR
  #define PROXY

  #include "../../config.h"
  #include "../../misc/macros.hpp"
  #include "../../traits/realTraits.hpp"
  #include "../expressionInterface.hpp"
  #define OPERATOR ==
  #define PASSIVE_TYPE double

namespace codi {
#endif

  /// Function overload for operator OPERATOR.
  template<typename Real, typename ArgA>
  CODI_INLINE bool operator OPERATOR(ExpressionInterface<Real, ArgA> const& argA, PASSIVE_TYPE const& argB) {
    return RealTraits::getPassiveValue(argA.cast()) OPERATOR argB;
  }

  /// Function overload for operator OPERATOR.
  template<typename Real, typename ArgB>
  CODI_INLINE bool operator OPERATOR(PASSIVE_TYPE const& argA, ExpressionInterface<Real, ArgB> const& argB) {
    return argA OPERATOR RealTraits::getPassiveValue(argB.cast());
  }

// Create a correct include environment for viewing and programming in an IDE.
#ifdef PROXY
  #undef PROXY
}
#endif

#undef PASSIVE_TYPE
