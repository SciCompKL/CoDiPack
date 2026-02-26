/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2026 Chair for Scientific Computing (SciComp), RPTU University Kaiserslautern-Landau
 * Homepage: http://scicomp.rptu.de
 * Contact:  Prof. Nicolas R. Gauger (codi@scicomp.uni-kl.de)
 *
 * Lead developers: Max Sagebaum, Johannes Blühdorn (SciComp, RPTU University Kaiserslautern-Landau)
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
 *  - SciComp, RPTU University Kaiserslautern-Landau:
 *    - Max Sagebaum
 *    - Johannes Blühdorn
 *    - Former members:
 *      - Tim Albring
 */
#pragma once

#include "expressionTraits.hpp"
#include "../expressions/activeType.hpp"

#include <Eigen/src/Core/MathFunctions.h>

namespace Eigen {
  namespace internal {

    /// Required for eigen math implementations that are called with expressions of CoDiPack types, e.g., sqrt(a+b);
    /// The definition is not optimal since we cut of the expression. If want to do this in an optimal way, we would
    /// have to specialize all math implementation functions, e.g. Eigen::internal::sqrt_impl.
    template<typename T>
    struct global_math_functions_filtering_base<T, codi::ExpressionTraits::EnableIfExpression<T>> {

        using type = codi::ActiveType<typename T::ADLogic>; ///< Not optimal since we cut of the expression.
    };
  }
}
