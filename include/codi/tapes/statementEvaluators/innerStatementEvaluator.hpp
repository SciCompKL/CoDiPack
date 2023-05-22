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

#include <algorithm>
#include <functional>
#include <type_traits>

#include "../../expressions/activeType.hpp"
#include "../../misc/macros.hpp"
#include "../../traits/expressionTraits.hpp"
#include "directStatementEvaluator.hpp"
#include "statementEvaluatorInterface.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Additional data required by an InnerStatementEvaluator.
   */
  struct InnerPrimalTapeStatementData : public PrimalTapeStatementFunctions {
    public:

      using Base = PrimalTapeStatementFunctions;  ///< Base class abbreviation.

      size_t maxActiveArguments;    ///< Maximum number of active arguments.
      size_t maxConstantArguments;  ///< Maximum number of constant arguments.

      /// Constructor
      InnerPrimalTapeStatementData(size_t maxActiveArguments, size_t maxConstantArguments,
                                   typename Base::Handle forward, typename Base::Handle primal,
                                   typename Base::Handle reverse)
          : Base(forward, primal, reverse),
            maxActiveArguments(maxActiveArguments),
            maxConstantArguments(maxConstantArguments) {}
  };

  /// Store InnerPrimalTapeStatementData as static variables for each combination of generator (tape) and expression
  /// used in the program.
  template<typename Tape, typename Expr>
  struct InnerStatementEvaluatorStaticStore {
    public:

      /// Static storage. Static construction is done by instantiating the statementEvaluate*Inner functions of the
      /// generator with Expr. Also evaluates the number of active type arguments and constant type arguments.
      static InnerPrimalTapeStatementData const staticStore;  ///< Static storage.
  };

  template<typename Generator, typename Expr>
  InnerPrimalTapeStatementData const InnerStatementEvaluatorStaticStore<Generator, Expr>::staticStore(
      ExpressionTraits::NumberOfActiveTypeArguments<Expr>::value,
      ExpressionTraits::NumberOfConstantTypeArguments<Expr>::value,
      (typename PrimalTapeStatementFunctions::Handle)Generator::template statementEvaluateForwardInner<Expr>,
      (typename PrimalTapeStatementFunctions::Handle)Generator::template statementEvaluatePrimalInner<Expr>,
      (typename PrimalTapeStatementFunctions::Handle)Generator::template statementEvaluateReverseInner<Expr>);

  /**
   * @brief Expression evaluation in the inner function. Data loading in the compilation context of the tape.
   * Storing in static context.
   *
   * Data loading is performed in the compilation context of the tape. The tape will then call the handle for the
   * evaluation of the expression after the data is loaded. This evaluator stores expression specific data and the
   * inner function handles.
   *
   * See StatementEvaluatorInterface for details.
   *
   * @tparam T_Real  The computation type of a tape, usually chosen as ActiveType::Real.
   */
  template<typename T_Real>
  struct InnerStatementEvaluator : public StatementEvaluatorInterface<T_Real> {
    public:

      using Real = CODI_DD(T_Real, double);  ///< See InnerStatementEvaluator.

      /*******************************************************************************/
      /// @name StatementEvaluatorInterface implementation
      /// @{

      using Handle = InnerPrimalTapeStatementData const*;  ///< Pointer to static storage location.

      /// \copydoc StatementEvaluatorInterface::callForward
      template<typename Tape, typename... Args>
      static Real callForward(Handle const& h, Args&&... args) {
        return Tape::statementEvaluateForwardFull((FunctionForward<Tape>)h->forward, h->maxActiveArguments,
                                                  h->maxConstantArguments, std::forward<Args>(args)...);
      }

      /// \copydoc StatementEvaluatorInterface::callPrimal
      template<typename Tape, typename... Args>
      static Real callPrimal(Handle const& h, Args&&... args) {
        return Tape::statementEvaluatePrimalFull((FunctionPrimal<Tape>)h->primal, h->maxActiveArguments,
                                                 h->maxConstantArguments, std::forward<Args>(args)...);
      }

      /// \copydoc StatementEvaluatorInterface::callReverse
      template<typename Tape, typename... Args>
      static void callReverse(Handle const& h, Args&&... args) {
        Tape::statementEvaluateReverseFull((FunctionReverse<Tape>)h->reverse, h->maxActiveArguments,
                                           h->maxConstantArguments, std::forward<Args>(args)...);
      }

      /// \copydoc StatementEvaluatorInterface::createHandle
      template<typename Tape, typename Generator, typename Expr>
      static Handle createHandle() {
        return &InnerStatementEvaluatorStaticStore<Generator, Expr>::staticStore;
      }

      /// @}

    protected:

      /// Full forward function type.
      template<typename Tape>
      using FunctionForward = decltype(&Tape::template statementEvaluateForwardInner<ActiveType<Tape>>);

      /// Full primal function type.
      template<typename Tape>
      using FunctionPrimal = decltype(&Tape::template statementEvaluatePrimalInner<ActiveType<Tape>>);

      /// Full reverse function type.
      template<typename Tape>
      using FunctionReverse = decltype(&Tape::template statementEvaluateReverseInner<ActiveType<Tape>>);
  };
}
