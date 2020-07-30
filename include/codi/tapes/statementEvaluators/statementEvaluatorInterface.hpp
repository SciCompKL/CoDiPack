#pragma once

#include "../../aux/macros.hpp"
#include "../../aux/memberStore.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  template<typename _Real>
  struct StatementEvaluatorInterface {

      using Real = CODI_DECLARE_DEFAULT(_Real, double);

    public:
      /*******************************************************************************
       * Section: Start of interface definition
       *
       */
      using Handle = CODI_ANY;

      template<typename Tape, typename ... Args>
      static Real callForward(Handle const& h, Args&& ... args);

      template<typename Tape, typename ... Args>
      static Real callPrimal(Handle const& h, Args&& ... args);

      template<typename Tape, typename ... Args>
      static void callReverse(Handle const& h, Args&& ... args);

      template<typename Tape, typename Generator, typename Expr>
      static Handle createHandle();

  };
}
