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

#include <algorithm>
#include <functional>
#include <type_traits>

#include "../../misc/macros.hpp"
#include "../../expressions/activeType.hpp"
#include "../../expressions/assignExpression.hpp"
#include "../../traits/expressionTraits.hpp"
#include "../misc/statementSizes.hpp"
#include "directStatementEvaluator.hpp"
#include "statementEvaluatorInterface.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Additional data required by an InnerStatementEvaluator.
   */
  struct InnerPrimalTapeStatementData {
    public:

      using Base = PrimalTapeStatementFunctions;  ///< Base class abbreviation.

      PrimalTapeStatementFunctions functions;
      StatementSizes stmtSizes;

      /// Constructor
      InnerPrimalTapeStatementData(PrimalTapeStatementFunctions functions, StatementSizes stmtSizes)
          : functions(functions),
            stmtSizes(stmtSizes) {}
  };

  /// Store InnerPrimalTapeStatementData as static variables for each combination of generator (tape) and expression
  /// used in the program.
  template<typename Generator, typename Expr>
  struct InnerStatementEvaluatorStaticStore {
    public:

      /// Static storage. Static construction is done by instantiating the statementEvaluate*Inner functions of the
      /// generator with Expr. Also evaluates the number of active type arguments and constant type arguments.
      static InnerPrimalTapeStatementData const staticStore;  ///< Static storage.

      /// Generates the data for the static store.
      template<StatementCall ... types>
      static InnerPrimalTapeStatementData gen() {
        using Handle = typename PrimalTapeStatementFunctions::Handle;

        return InnerPrimalTapeStatementData(
            PrimalTapeStatementFunctions(
              ((Handle)Generator::template StatementCallGen<types, Expr>::evaluateInner)...
            ),
            StatementSizes::create<Expr>());
      }
  };

  template<typename Generator, typename Expr>
  InnerPrimalTapeStatementData const InnerStatementEvaluatorStaticStore<Generator, Expr>::staticStore =
      InnerStatementEvaluatorStaticStore<Generator, Expr>::gen<CODI_STMT_CALL_GEN_ARGS>();

  /**
   * @brief Expression evaluation in the inner function. Data loading in the compilation context of the tape.
   * Storing in static context.
   *
   * Data loading is performed in the compilation context of the tape. The tape will then call the handle for the
   * evaluation of the expression after the data is loaded. This evaluator stores expression specific data and the
   * inner function handles.
   *
   * See StatementEvaluatorInterface for details.
   */
  struct InnerStatementEvaluator : public StatementEvaluatorInterface {
    public:

      /*******************************************************************************/
      /// @name StatementEvaluatorInterface implementation
      /// @{

      using Handle = InnerPrimalTapeStatementData const*;  ///< Pointer to static storage location.

      /// \copydoc StatementEvaluatorInterface::call
      template<StatementCall type, typename Tape, typename... Args>
      static void call(Handle const& h, Args&&... args) {
        using Expr = AssignExpression<ActiveType<Tape>, ActiveType<Tape>>;
        using CallGen = typename Tape::template StatementCallGen<type, Expr>;

        using Function = decltype(&CallGen::evaluateInner);

        CallGen::evaluateFull(((Function)h->functions.funcs[(size_t)type]), h->stmtSizes, std::forward<Args>(args)...);
      }

      /// \copydoc StatementEvaluatorInterface::createHandle
      template<typename Tape, typename Generator, typename Expr>
      static Handle createHandle() {
        return &InnerStatementEvaluatorStaticStore<Generator, Expr>::staticStore;
      }

      /// @}
  };
}
