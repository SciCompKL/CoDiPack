#pragma once

#include "../../aux/macros.h"
#include "../../aux/memberStore.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  template<typename _Real>
  struct StatementEvaluatorInterface {

      using Real = DECLARE_DEFAULT(_Real, double);

    public:
      /*******************************************************************************
       * Section: Start of interface definition
       *
       */
      using Handle = ANY;

      template<typename Tape, typename Expr>
      static Handle createHandle();

      template<typename Tape, typename ... Args>
      static void callReverse(Handle const& h, Args&& ... args);

      template<typename Tape, typename ... Args>
      static Real callForward(Handle const& h, Args&& ... args);

      template<typename Tape, typename ... Args>
      static Real callPrimal(Handle const& h, Args&& ... args);
  };
}
