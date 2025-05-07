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

#include "../../misc/macros.hpp"
#include "../../misc/memberStore.hpp"
#include "../misc/statementSizes.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /// Defines all the operations which can be evaluated on a statement by a tape.
  enum class StatementCall {
    ClearAdjoints,     ///< Clear the adjoint values.
    Forward,           ///< Evaluate expression in a forward mode.
    Primal,            ///< Evaluate primal expression.
    ResetPrimals,      ///< Restore the primal values.
    Reverse,           ///< Evaluate expression in a reverse mode.
    WriteInformation,  ///< Get write information.
    N_Elements         ///< Number of elements.
  };

#define CODI_STMT_CALL_GEN_ARGS                                                                             \
  StatementCall::ClearAdjoints, StatementCall::Forward, StatementCall::Primal, StatementCall::ResetPrimals, \
      StatementCall::Reverse, StatementCall::WriteInformation

  /**
   * @brief Tape side interface for StatementEvaluatorInterface.
   *
   * See StatementEvaluatorInterface for a full description.
   *
   * In the evaluation method of StatementCallGenerator the full evaluation of the statement needs to be done.
   * - 1. Load expression specific data
   * - 2. Call expression specific function
   *
   * The structure StatementCallGenerator needs to be specialized for all enums in StatementCall.
   */
  struct StatementEvaluatorTapeInterface {
    public:

      /*******************************************************************************/
      /// @name Interface definition

      /// This structure is accessed by the StatementEvaluatorInterface.
      template<StatementCall type, typename Expr>
      struct StatementCallGenerator {
          /// Evaluate the full expression.
          template<typename... Args>
          CODI_INLINE static void evaluate(Args&&... args);
      };
  };

  /**
   * @brief Tape side interface for StatementEvaluatorInterface.
   *
   * See StatementEvaluatorInterface for a full description.
   *
   * The `evaluateInner` methods need to be stored by the StatementEvaluatorInterface. These methods
   * perform the `Call expression specific function` logic.
   *
   * The `evaluateFull` functions are called by the StatementEvaluatorInterface on a `call` function call.
   * This performs the step `Load expression specific data` in an inline context. `inner` is the stored function pointer
   * in the handle.
   */
  struct StatementEvaluatorInnerTapeInterface {
    public:

      /*******************************************************************************/
      /// @name Interface definition

      /// This structure is accessed by the StatementEvaluatorInterface.
      template<StatementCall type, typename Expr>
      struct StatementCallGenerator {
          /// Evaluate expression in a forward mode.
          template<typename... Args>
          CODI_INLINE static void evaluateInner(Args&&... args);

          /// Load the expression data and evaluate the expression.
          template<typename InnerFunc, typename... Args>
          CODI_INLINE static void evaluateFull(InnerFunc func, Args&&... args);
      };
  };
}
