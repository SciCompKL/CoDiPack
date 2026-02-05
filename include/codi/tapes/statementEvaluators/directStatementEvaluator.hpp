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

#include <algorithm>
#include <functional>
#include <type_traits>

#include "../../expressions/activeType.hpp"
#include "../../misc/macros.hpp"
#include "../misc/assignStatement.hpp"
#include "statementEvaluatorInterface.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Data required for all possible handle calls.
   */
  struct PrimalTapeStatementFunctions {
    public:

      using Handle = void*;  ///< Function pointer.

      std::array<Handle, (size_t)StatementCall::N_Elements> funcs;  ///< Array for the function handles.

      /// Constructor
      template<typename... Args>
      PrimalTapeStatementFunctions(Args... args) : funcs{args...} {}
  };

  /// Store PrimalTapeStatementFunctions as static variables for each combination of generator (tape) and expression
  /// used in the program.
  template<typename Generator, typename Stmt>
  struct DirectStatementEvaluatorStaticStore {
    public:

      /// Static storage. Static construction is done by instantiating the statementEvaluate* functions of the generator
      /// with Stmt.
      static PrimalTapeStatementFunctions const staticStore;

      /// Generates the data for the static store.
      template<StatementCall... types>
      static PrimalTapeStatementFunctions gen() {
        using Handle = typename PrimalTapeStatementFunctions::Handle;

        return PrimalTapeStatementFunctions(
            reinterpret_cast<Handle>(Generator::template StatementCallGenerator<types, Stmt>::evaluate)...);
      }
  };

  template<typename Generator, typename Stmt>
  PrimalTapeStatementFunctions const DirectStatementEvaluatorStaticStore<Generator, Stmt>::staticStore =
      DirectStatementEvaluatorStaticStore<Generator, Stmt>::gen<CODI_STMT_CALL_GEN_ARGS>();

  /**
   * @brief Full evaluation of the expression in the function handle. Storing in static context.
   *
   * Data loading and evaluation of the expression are all done in the handle. This evaluator will directly evaluate the
   * full handle for the expression.
   *
   * See StatementEvaluatorInterface for details.
   */
  struct DirectStatementEvaluator : public StatementEvaluatorInterface {
    public:

      /*******************************************************************************/
      /// @name StatementEvaluatorInterface implementation
      /// @{

      using Handle = PrimalTapeStatementFunctions const*;  ///< Pointer to static storage location.

      /// \copydoc StatementEvaluatorInterface::call
      template<StatementCall type, typename Tape, typename... Args>
      static void call(Handle const& h, Args&&... args) {
        using Stmt = AssignStatement<ActiveType<Tape>, ActiveType<Tape>>;
        using CallGen = typename Tape::template StatementCallGenerator<type, Stmt>;

        using Function = decltype(&CallGen::evaluate);

        ((Function)h->funcs[(size_t)type])(std::forward<Args>(args)...);
      }

      /// \copydoc StatementEvaluatorInterface::createHandle
      template<typename Tape, typename Generator, typename Stmt>
      static Handle createHandle() {
        return &DirectStatementEvaluatorStaticStore<Generator, Stmt>::staticStore;
      }

      /// @}
  };
}
