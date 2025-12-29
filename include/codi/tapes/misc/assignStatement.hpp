/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2026 Chair for Scientific Computing (SciComp), University of Kaiserslautern-Landau
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
#include "../../expressions/lhsExpressionInterface.hpp"
#include "../../misc/macros.hpp"

/** \copydoc codi::Namespace */
namespace codi {
  /**
   * Pseudo expression that defines an lhs and rhs expression where the rhs is assigned to the lhs.
   *
   * @tparam T_Lhs The type of the left hand side expression that is assigned to.
   * @tparam T_Rhs The type of the right hand side expression that is assigned from.
   */
  template<typename T_Lhs, typename T_Rhs>
  struct AssignStatement {
      using Lhs =
          CODI_DD(T_Lhs, CODI_T(LhsExpressionInterface<double, double, CODI_ANY, CODI_ANY>));  ///< See AssignStatement.
      using Rhs = CODI_DD(T_Rhs, CODI_T(ExpressionInterface<double, CODI_ANY>));               ///< See AssignStatement.
  };
}
