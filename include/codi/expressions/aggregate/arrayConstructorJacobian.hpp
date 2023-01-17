/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2022 Chair for Scientific Computing (SciComp), TU Kaiserslautern
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

#include "../../config.h"
#include "../../misc/macros.hpp"
#include "../../traits/realTraits.hpp"
#include "../expressionInterface.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /// Helper class for the detection of a Jacobian application for an array constructor.
  /// @tparam T_Creator The expression creating the Jacobian.
  /// @tparam T_ReturnType Used in the multiplication detection for the return type.
  template<typename T_Creator, typename T_ReturnType, size_t T_index>
  struct ArrayConstructorJacobian {
      using Creator = CODI_DD(T_Creator, CODI_T(ExpressionInterface<double, void>));  ///< See ArrayConstructorJacobian.
      using ReturnType = CODI_DD(T_ReturnType, double);                               ///< See ArrayConstructorJacobian.

      static size_t constexpr index = CODI_DD(T_index, 0);  ///< The index that is accessed.

      Creator const& creator;  ///< Reference to the creator.

      /// Constructor.
      CODI_INLINE ArrayConstructorJacobian(Creator const& creator) : creator(creator) {}
  };

  /// Detection of the application of a Jacobian from an array constructor. See ArrayConstructorJacobian.
  template<typename Type, typename Creator, typename ReturnType, size_t index>
  CODI_INLINE ReturnType operator*(ArrayConstructorJacobian<Creator, ReturnType, index> const& reduce,
                                   Type const& jac) {
    return RealTraits::AggregatedTypeTraits<typename Creator::Real>::template adjointOfConstructor<index>(
        reduce.creator.getValue(), jac);
  }
}
