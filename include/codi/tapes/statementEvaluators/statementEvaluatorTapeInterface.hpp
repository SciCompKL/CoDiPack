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
#pragma once

#include "../../misc/macros.hpp"
#include "../../misc/memberStore.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Tape side interface for StatementEvaluatorInterface.
   *
   * See StatementEvaluatorInterface for a full description.
   *
   * In every method the full evaluation of the statement needs to be done.
   * - 1. Load expression specific data
   * - 2. Call expression specific function
   *
   * @tparam T_Real  The computation type of a tape usually defined by ActiveType::Real.
   */
  template<typename T_Real>
  struct StatementEvaluatorTapeInterface {
    public:

      using Real = CODI_DD(T_Real, double);  ///< See StatementEvaluatorTapeInterface

      /*******************************************************************************/
      /// @name Interface definition

      /// Evaluate expression in a forward mode.
      template<typename Expr, typename... Args>
      static Real statementEvaluateForward(Args&&... args);

      /// Evaluate primal expression.
      template<typename Expr, typename... Args>
      static Real statementEvaluatePrimal(Args&&... args);

      /// Evaluate expression in a reverse mode.
      template<typename Expr, typename... Args>
      static void statementEvaluateReverse(Args&&... args);
  };

  /**
   * @brief Tape side interface for StatementEvaluatorInterface.
   *
   * See StatementEvaluatorInterface for a full description.
   *
   * The `statementEvaluate*Inner` methods needs to be stored by the StatementEvaluatorInterface. These methods
   * perform the `Call expression specific function` logic.
   *
   * The `statementEvaluate*Full` functions are called by the StatementEvaluatorInterface on a `call*` function call.
   * This performs the step `Load expression specific data` in an inline context. `inner` is the stored function pointer
   * in the handle.
   *
   * @tparam T_Real  The computation type of a tape usually defined by ActiveType::Real.
   */
  template<typename T_Real>
  struct StatementEvaluatorInnerTapeInterface {
    public:

      using Real = CODI_DD(T_Real, double);  ///< See StatementEvaluatorInnerTapeInterface

      /*******************************************************************************/
      /// @name Interface definition

      /// Load the expression data and evaluate the expression in a forward mode.
      template<typename Func, typename... Args>
      static Real statementEvaluateForwardFull(Func const& inner, size_t const& maxActiveArgs,
                                               size_t const& maxConstantArgs, Args&&... args);

      /// Load the expression data and evaluate the expression in a primal setting.
      template<typename Func, typename... Args>
      static Real statementEvaluatePrimalFull(Func const& inner, size_t const& maxActiveArgs,
                                              size_t const& maxConstantArgs, Args&&... args);

      /// Load the expression data and evaluate the expression in a reverse mode.
      template<typename Func, typename... Args>
      static void statementEvaluateReverseFull(Func const& inner, size_t const& maxActiveArgs,
                                               size_t const& maxConstantArgs, Args&&... args);

      /// Evaluate expression in a forward mode.
      template<typename Expr, typename... Args>
      static Real statementEvaluateForwardInner(Args&&... args);

      /// Evaluate expression in a primal setting.
      template<typename Expr, typename... Args>
      static Real statementEvaluatePrimalInner(Args&&... args);

      /// Evaluate expression in a reverse mode.
      template<typename Expr, typename... Args>
      static void statementEvaluateReverseInner(Args&&... args);
  };
}
