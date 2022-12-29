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
   * @brief Creation of handles for the evaluation of expressions in a context where the expression type is not
   * available.
   *
   * For primal value taping, the expression templates create a problem. The expression type and value is generated in
   * the code where the statement is evaluated. The type and value of the expression is required again to evaluate
   * the reverse or forward mode of the elemental function described by the expression (see
   * \ref AD_TheoryMathematicalDefinitions). The problem is that types cannot be stored in C++. Also, the value of the
   * expression is context specific and this context is not available during the tape interpretation either. The
   * generation of handles solves this problem by providing the tape with functions for the evaluation of the
   * expressions. Since these functions can be instantiated with the type of the expression, inside of the handles, the
   * expression type is available again. For more details on the topic please see \ref SAG2018Expression.
   *
   * The creation of the handles is tightly coupled with the implementing tape. The tape defines how the data for the
   * expression is loaded and which functions are evaluated on the expression in order to perform the necessary
   * operations for the tape. However, how the handles are generated, how the data for the handles is stored (e.g. which
   * kinds of function pointers) and how they are evaluated is neither tape specific nor is one way optimal for every
   * use case.
   *
   * The process of evaluating a handle can be separated into multiple steps:
   * - 1. Load statement specific data
   * - 2. Load expression specific data
   * - 3. Call expression specific function
   *
   * The call to a handle may take place either between step 1 and 2 or 2 and 3. The advantage between step 1 and 2 is
   * that it is very simple to implement, the disadvantage is that the expression specific data load is handled after
   * the function pointer call. The compiler can therefore not optimize these data loading calls. If the call of the
   * handle is between step 2 and 3, then the compiler can optimize the generalized data handling and only the specifics
   * for the expression are hidden behind the function pointer call. However, in order to perform step 2, the tape needs
   * some expression specific data which has to be extracted from the handle. Therefore, the implementation is more
   * involved.
   *
   * The first approach (handle call between step 1 and 2) is defined by the StatementEvaluatorTapeInterface.
   * The second approach (handle call between step 2 and 3) is defined by the StatementEvaluatorInnerTapeInterface.
   *
   * In general, the tape implements the interfaces mentioned above. However, special cases like preaccumulation support
   * requires a different generator. Therefore, `createHandle` has two template arguments, one for the tape and one for
   * the generator.
   *
   * In general, implementations of this interface need to store functions pointers to the `statementEvaluate*`
   * functions of the StatementEvaluatorTapeInterface or function pointers to the `statementEvaluate*Inner` of the
   * StatementEvaluatorInnerTapeInterface.
   *
   * A usual call flow for the first approach is (see also the code for ReverseStatementEvaluator):
   * \code{.cpp}
   *   // During recording in tape.store.
   *   // Instantiates e.g Tape::statementEvaluateReverse<Expr>
   *   auto handle = StatementEvaluatorInterface::createHandle<Tape, Tape, Expr>();
   *   tapeData.pushHandle(handle);
   *
   *   // During reverse interpretation of the tape.
   *   auto handle = tapeData.popHandle(handle);
   *
   *   // This calls Tape::statementEvaluateReverse<Expr>(tapeData);
   *   StatementEvaluatorInterface::callReverse(handle, tapeData);
   * \endcode
   *
   * @tparam T_Real  The computation type of a tape, usually chosen as ActiveType::Real.
   */
  template<typename T_Real>
  struct StatementEvaluatorInterface {
    public:

      using Real = CODI_DD(T_Real, double);  ///< See StatementEvaluatorInterface.

      /*******************************************************************************/
      /// @name Interface definition

      using Handle = CODI_ANY;  ///< Type of the handle.

      /// @tparam Tape  Has to implement StatementEvaluatorTapeInterface or StatementEvaluatorInnerTapeInterface,
      ///               depending on the interface the implementation uses.
      template<typename Tape, typename... Args>
      static Real callForward(Handle const& h, Args&&... args);

      /// @tparam Tape  Has to implement StatementEvaluatorTapeInterface or StatementEvaluatorInnerTapeInterface,
      ///               depending on the interface the implementation uses.
      template<typename Tape, typename... Args>
      static Real callPrimal(Handle const& h, Args&&... args);

      /// @tparam Tape  Has to implement StatementEvaluatorTapeInterface or StatementEvaluatorInnerTapeInterface,
      ///               depending on the interface the implementation uses.
      template<typename Tape, typename... Args>
      static void callReverse(Handle const& h, Args&&... args);

      /// @tparam Tape       Usually not required. Access tape specific configurations.
      /// @tparam Generator  Has to implement the StatementEvaluatorTapeInterface or
      ///                    StatementEvaluatorInnerTapeInterface. Usually the same as Tape.
      /// @tparam Expr       Instance of ExpressionInterface.
      template<typename Tape, typename Generator, typename Expr>
      static Handle createHandle();
  };
}
