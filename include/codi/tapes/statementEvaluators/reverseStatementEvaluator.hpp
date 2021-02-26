#pragma once

#include <algorithm>
#include <functional>
#include <type_traits>

#include "../../aux/macros.hpp"
#include "../../aux/memberStore.hpp"
#include "../../aux/exceptions.hpp"
#include "../../expressions/activeType.hpp"
#include "statementEvaluatorInterface.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Only stored the function handle for the reverse evaluation.
   *
   * Uses the StatementEvaluatorTapeInterface.
   *
   * @tparam _Real  The computation type of a tape usually defined by ActiveType::Real.
   */
  template<typename _Real>
  struct ReverseStatementEvaluator : public StatementEvaluatorInterface<_Real> {

      using Real = CODI_DD(_Real, double);  ///< See ReverseStatementEvaluator

    public:

      /*******************************************************************************/
      /// @name StatementEvaluatorInterface implementation
      /// @{

      using Handle = void*;  ///< Function pointer to the reverse evaluation

      /// Throws CODI_EXCEPTION on call.
      template<typename Tape, typename... Args>
      static Real callForward(Handle const& h, Args&&... args) {
        CODI_UNUSED(h, args...);

        CODI_EXCEPTION("ReverseStatementEvaluator does not support forward evaluation calls.");

        return Real();
      }

      /// Throws CODI_EXCEPTION on call.
      template<typename Tape, typename... Args>
      static Real callPrimal(Handle const& h, Args&&... args) {
        CODI_UNUSED(h, args...);

        CODI_EXCEPTION("ReverseStatementEvaluator does not support primal evaluation calls.");

        return Real();
      }

      /// \copydoc StatementEvaluatorInterface::callReverse
      template<typename Tape, typename... Args>
      static void callReverse(Handle const& h, Args&&... args) {
        HandleTyped<Tape> func = (HandleTyped<Tape>)h;

        func(std::forward<Args>(args)...);
      }

      /// \copydoc StatementEvaluatorInterface::createHandle
      template<typename Tape, typename Generator, typename Expr>
      static Handle createHandle() {
        return (Handle*)Generator::template statementEvaluateReverse<Expr>;
      }

      /// @}

    protected:

      /// Full reverse function type.
      template<typename Tape>
      using HandleTyped = decltype(&Tape::template statementEvaluateReverse<ActiveType<Tape>>);

  };
}
