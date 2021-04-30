#pragma once

#include <algorithm>
#include <functional>
#include <type_traits>

#include "../../aux/macros.hpp"
#include "../../expressions/activeType.hpp"
#include "statementEvaluatorInterface.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Data required for all possible handle calls.
   */
  struct PrimalTapeStatementFunctions {
    public:

      using Handle = void*;  ///< Function pointer.

      Handle forward;  ///< Forward function handle.
      Handle primal;   ///< Primal function handle.
      Handle reverse;  ///< Reverse function handle.

      /// Constructor
      PrimalTapeStatementFunctions(Handle forward, Handle primal, Handle reverse)
          : forward(forward), primal(primal), reverse(reverse) {}
  };

  /// Store PrimalTapeStatementFunctions as static variables for each combination of generator (tape) and expression
  /// used in the program.
  template<typename Generator, typename Expr>
  struct DirectStatementEvaluatorStaticStore {
    public:

      /// Static storage. Static construction is done by instantiating the statementEvaluate* functions of the generator
      /// with Expr.
      static PrimalTapeStatementFunctions const staticStore;
  };

  template<typename Generator, typename Expr>
  PrimalTapeStatementFunctions const DirectStatementEvaluatorStaticStore<Generator, Expr>::staticStore(
      (typename PrimalTapeStatementFunctions::Handle)Generator::template statementEvaluateForward<Expr>,
      (typename PrimalTapeStatementFunctions::Handle)Generator::template statementEvaluatePrimal<Expr>,
      (typename PrimalTapeStatementFunctions::Handle)Generator::template statementEvaluateReverse<Expr>);

  /**
   * @brief Full evaluation of the expression in the function handle. Storing in static context.
   *
   * Data loading and evaluation of the expression are all done in the handle. This evaluator will directly evaluate the
   * full handle for the expression.
   *
   * See StatementEvaluatorInterface for details.
   *
   * @tparam _Real  The computation type of a tape, usually chosen as ActiveType::Real.
   */
  template<typename _Real>
  struct DirectStatementEvaluator : public StatementEvaluatorInterface<_Real> {
    public:

      using Real = CODI_DD(_Real, double);  ///< See DirectStatementEvaluator.

      /*******************************************************************************/
      /// @name StatementEvaluatorInterface implementation
      /// @{

      using Handle = PrimalTapeStatementFunctions const*;  ///< Pointer to static storage location.

      /// \copydoc StatementEvaluatorInterface::callForward
      template<typename Tape, typename... Args>
      static Real callForward(Handle const& h, Args&&... args) {
        return ((FunctionForward<Tape>)h->forward)(std::forward<Args>(args)...);
      }

      /// \copydoc StatementEvaluatorInterface::callPrimal
      template<typename Tape, typename... Args>
      static Real callPrimal(Handle const& h, Args&&... args) {
        return ((FunctionPrimal<Tape>)h->primal)(std::forward<Args>(args)...);
      }

      /// \copydoc StatementEvaluatorInterface::callReverse
      template<typename Tape, typename... Args>
      static void callReverse(Handle const& h, Args&&... args) {
        ((FunctionReverse<Tape>)h->reverse)(std::forward<Args>(args)...);
      }

      /// \copydoc StatementEvaluatorInterface::createHandle
      template<typename Tape, typename Generator, typename Expr>
      static Handle createHandle() {
        return &DirectStatementEvaluatorStaticStore<Generator, Expr>::staticStore;
      }

      /// @}

    protected:

      /// Full forward function type.
      template<typename Tape>
      using FunctionForward = decltype(&Tape::template statementEvaluateForward<ActiveType<Tape>>);

      /// Full primal function type.
      template<typename Tape>
      using FunctionPrimal = decltype(&Tape::template statementEvaluatePrimal<ActiveType<Tape>>);

      /// Full reverse function type.
      template<typename Tape>
      using FunctionReverse = decltype(&Tape::template statementEvaluateReverse<ActiveType<Tape>>);
  };
}
