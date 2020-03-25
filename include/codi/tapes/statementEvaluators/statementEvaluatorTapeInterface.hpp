#pragma once

#include "../../aux/macros.h"
#include "../../aux/memberStore.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  template<typename _Real>
  struct StatementEvaluatorTapeInterface {
    public:

      using Real = DECLARE_DEFAULT(_Real, double);

      /*******************************************************************************
       * Section: Start of interface definition
       *
       */

      template<typename Expr, typename ... Args>
      static void statementEvaluateReverse(Args&& ... args);

      template<typename Expr, typename ... Args>
      static Real statementEvaluateForward(Args&& ... args);

      template<typename Expr, typename ... Args>
      static Real statementEvaluatePrimal(Args&& ... args);
  };
}
