#pragma once

#include "../../aux/macros.hpp"
#include "../../aux/memberStore.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Creation of handles for the evaluation of expressions in a different context. (Where the expression type
   * is not available.)
   *
   * The creation of the handles is tightly coupled with the implementing tape. The tape defines how the data for the
   * expression is loaded and which functions are evaluated on the expression in order to perform the necessary
   * operations for the tape. But how the handles are generated, how the data for the handles is stored and how they are
   * evaluated is neither tape specific nor is one way optimal for every use case.
   *
   * The process of evaluating a handle can be separated into multiple steps:
   * - 1. Load statement specific data
   * - 2. Load expression specific data
   * - 3. Call expression specific function
   *
   * The call to a handle can now either be between step 1 and 2 or 2 and 3. The advantage between step 1. and 2. is that
   * it is very simple to implement, the disadvantage is that the expression specific data load is handled after the function
   * pointer call. The compiler can therefore not optimize these data loading calls. If the call of the handle is
   * between step 2 and 3, then the compiler can optimize the generalized data handling and only the specifics for the
   * expression are hidden behind the function pointer call. But in order to perform step 2 the tape needs some
   * expression specific data which need to be extracted from the handle. Therefore the implementation is more involved.
   *
   * The first approach, handle call between step 1 and 2 is defined by the StatementEvaluatorTapeInterface.
   * The second approach, handle call between step 2 and 3 is defined by the StatementEvaluatorInnerTapeInterface.
   *
   * In general, implementations of this interface need to store functions pointers to the `statementEvaluate*` functions
   * of the StatementEvaluatorTapeInterface or function pointers to the `statementEvaluateForwardInner` of the
   * StatementEvaluatorInnerTapeInterface.
   *
   * @tparam _Real  The computation type of a tape usually defined by ActiveType::Real.
   */
  template<typename _Real>
  struct StatementEvaluatorInterface {
    public:

      using Real = CODI_DD(_Real, double); ///< See StatementEvaluatorInterface

      /*******************************************************************************/
      /// @name Interface definition

      using Handle = CODI_ANY;  ///< Type of the handle

      /// @tparam Tape  Needs to implement StatementEvaluatorTapeInterface or StatementEvaluatorInnerTapeInterface
      ///               depending which interface the implementation uses.
      template<typename Tape, typename ... Args>
      static Real callForward(Handle const& h, Args&& ... args);

      /// @tparam Tape  Needs to implement StatementEvaluatorTapeInterface or StatementEvaluatorInnerTapeInterface
      ///               depending which interface the implementation uses.
      template<typename Tape, typename ... Args>
      static Real callPrimal(Handle const& h, Args&& ... args);

      /// @tparam Tape  Needs to implement StatementEvaluatorTapeInterface or StatementEvaluatorInnerTapeInterface
      ///               depending which interface the implementation uses.
      template<typename Tape, typename ... Args>
      static void callReverse(Handle const& h, Args&& ... args);

      /// @tparam Tape       Usually not required. Access tape specific configurations.
      /// @tparam Generator  Needs to implement the StatementEvaluatorTapeInterface or StatementEvaluatorInnerTapeInterface
      /// @tparam Expr       Instance of ExpressionInterface.
      template<typename Tape, typename Generator, typename Expr>
      static Handle createHandle();

  };
}
