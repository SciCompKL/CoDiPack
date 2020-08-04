#pragma once

#include "../../aux/macros.hpp"
#include "../../aux/memberStore.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  template<typename _Real>
  struct StatementEvaluatorTapeInterface {
    public:

      using Real = CODI_DECLARE_DEFAULT(_Real, double);

      /*******************************************************************************
       * Section: Start of interface definition
       *
       */

      template<typename Expr, typename ... Args>
      static Real statementEvaluateForward(Args&& ... args);

      template<typename Expr, typename ... Args>
      static Real statementEvaluatePrimal(Args&& ... args);

      template<typename Expr, typename ... Args>
      static void statementEvaluateReverse(Args&& ... args);
  };

  template<typename _Real>
  struct StatementEvaluatorInnerTapeInterface {
    public:

      using Real = CODI_DECLARE_DEFAULT(_Real, double);

      /*******************************************************************************
       * Section: Start of interface definition
       *
       */

      template<typename Func, typename ... Args>
      static Real statementEvaluateForwardFull(
          Func const& inner, size_t const& maxActiveArgs, size_t const& maxConstantArgs, Args&& ... args);

      template<typename Func, typename ... Args>
      static Real statementEvaluatePrimalFull(
          Func const& inner, size_t const& maxActiveArgs, size_t const& maxConstantArgs, Args&& ... args);

      template<typename Func, typename ... Args>
      static void statementEvaluateReverseFull(
          Func const& inner, size_t const& maxActiveArgs, size_t const& maxConstantArgs, Args&& ... args);

      template<typename Expr, typename ... Args>
      static Real statementEvaluateForwardInner(Args&& ... args);

      template<typename Expr, typename ... Args>
      static Real statementEvaluatePrimalInner(Args&& ... args);

      template<typename Expr, typename ... Args>
      static void statementEvaluateReverseInner(Args&& ... args);
  };
}
