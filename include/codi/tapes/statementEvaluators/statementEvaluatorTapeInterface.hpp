#pragma once

#include "../../aux/macros.hpp"
#include "../../aux/memberStore.hpp"

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
   * @tparam _Real  The computation type of a tape usually defined by ActiveType::Real.
   */
  template<typename _Real>
  struct StatementEvaluatorTapeInterface {
    public:

      using Real = CODI_DECLARE_DEFAULT(_Real, double); ///< See StatementEvaluatorTapeInterface

      /*******************************************************************************/
      /// @name Interface definition

      /// Evaluate expression in a forward mode.
      template<typename Expr, typename ... Args>
      static Real statementEvaluateForward(Args&& ... args);

      /// Evaluate primal expression.
      template<typename Expr, typename ... Args>
      static Real statementEvaluatePrimal(Args&& ... args);

      /// Evaluate expression in a reverse mode.
      template<typename Expr, typename ... Args>
      static void statementEvaluateReverse(Args&& ... args);
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
   * @tparam _Real  The computation type of a tape usually defined by ActiveType::Real.
   */
  template<typename _Real>
  struct StatementEvaluatorInnerTapeInterface {
    public:

      using Real = CODI_DECLARE_DEFAULT(_Real, double); ///< See StatementEvaluatorInnerTapeInterface

      /*******************************************************************************/
      /// @name Interface definition

      /// Load the expression data and evaluate the expression in a forward mode.
      template<typename Func, typename ... Args>
      static Real statementEvaluateForwardFull(
          Func const& inner, size_t const& maxActiveArgs, size_t const& maxConstantArgs, Args&& ... args);

      /// Load the expression data and evaluate the expression in a primal setting.
      template<typename Func, typename ... Args>
      static Real statementEvaluatePrimalFull(
          Func const& inner, size_t const& maxActiveArgs, size_t const& maxConstantArgs, Args&& ... args);

      /// Load the expression data and evaluate the expression in a reverse mode.
      template<typename Func, typename ... Args>
      static void statementEvaluateReverseFull(
          Func const& inner, size_t const& maxActiveArgs, size_t const& maxConstantArgs, Args&& ... args);

      /// Evaluate expression in a forward mode.
      template<typename Expr, typename ... Args>
      static Real statementEvaluateForwardInner(Args&& ... args);

      /// Evaluate expression in a primal setting.
      template<typename Expr, typename ... Args>
      static Real statementEvaluatePrimalInner(Args&& ... args);

      /// Evaluate expression in a reverse mode.
      template<typename Expr, typename ... Args>
      static void statementEvaluateReverseInner(Args&& ... args);
  };
}
