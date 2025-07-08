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
#include "../../misc/macros.hpp"
#include "../../traits/expressionTraits.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /// Define all static sizes of an Expression.
  struct StatementSizes {
      size_t const outputArgs;    ///< Number of output arguments.
      size_t const inputArgs;     ///< Number of input arguments.
      size_t const constantArgs;  ///< Number of constant arguments.

      /// Constructor
      StatementSizes(size_t const outputArgs, size_t const inputArgs, size_t const constantArgs)
          : outputArgs(outputArgs), inputArgs(inputArgs), constantArgs(constantArgs) {}

      /// Creation function from AssignStatement.
      template<typename Stmt>
      static StatementSizes create() {
        using Rhs = typename Stmt::Rhs;
        using Lhs = typename Stmt::Lhs;

        size_t constexpr inputArgs = ExpressionTraits::NumberOfActiveTypeArguments<Rhs>::value;
        size_t constexpr constantArgs = ExpressionTraits::NumberOfConstantTypeArguments<Rhs>::value;
        size_t constexpr outputArgs = ExpressionTraits::NumberOfActiveTypeArguments<Lhs>::value;

        return StatementSizes(outputArgs, inputArgs, constantArgs);
      }
  };
}
