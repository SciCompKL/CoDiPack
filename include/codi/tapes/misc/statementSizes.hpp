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
#include "../../expressions/assignStatement.hpp"
#include "../../misc/macros.hpp"
#include "../../traits/expressionTraits.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /// Define all static sizes of an Expression.
  struct StatementSizes {
      size_t const inputArgs;     ///< Number of input arguments.
      size_t const constantArgs;  ///< Number of constant arguments.
      size_t const outputArgs;    ///< Number of output arguments.

      /// Constructor
      StatementSizes(size_t const inputArgs, size_t const constantArgs, size_t const maxOutputArgs)
          : inputArgs(inputArgs), constantArgs(constantArgs), outputArgs(maxOutputArgs) {}

      /// Creation function from AssignStatement.
      template<typename T_Stmt>
      static StatementSizes create() {
        using Stmt = CODI_DD(T_Stmt, CODI_T(AssignStatement<CODI_ANY, CODI_ANY>));

        using Lhs = typename Stmt::Lhs;
        using Rhs = typename Stmt::Rhs;

        size_t constexpr inputArgs = ExpressionTraits::NumberOfActiveTypeArguments<Rhs>::value;
        size_t constexpr constantArgs = ExpressionTraits::NumberOfConstantTypeArguments<Rhs>::value;
        size_t constexpr maxOutputArgs = ExpressionTraits::NumberOfActiveTypeArguments<Lhs>::value;

        return StatementSizes(inputArgs, constantArgs, maxOutputArgs);
      }
  };
}
