#pragma once

#include <algorithm>
#include <functional>
#include <type_traits>

#include "../../aux/macros.h"
#include "../../expressions/activeType.hpp"
#include "statementEvaluatorInterface.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  struct PrimalTapeStatementFunctions {
      using Handle = void*;

      Handle forward;
      Handle primal;
      Handle reverse;

      PrimalTapeStatementFunctions(Handle forward, Handle primal, Handle reverse) :
        forward(forward),
        primal(primal),
        reverse(reverse) {}
  };

  template<typename Generator, typename Expr>
  struct DirectStatementEvaluatorStaticStore {

      static PrimalTapeStatementFunctions const staticStore;
  };

  template<typename Generator, typename Expr>
  PrimalTapeStatementFunctions const DirectStatementEvaluatorStaticStore<Generator, Expr>::staticStore(
      (typename PrimalTapeStatementFunctions::Handle)Generator::template statementEvaluateForward<Expr>,
      (typename PrimalTapeStatementFunctions::Handle)Generator::template statementEvaluatePrimal<Expr>,
      (typename PrimalTapeStatementFunctions::Handle)Generator::template statementEvaluateReverse<Expr>);

  template<typename _Real>
  struct DirectStatementEvaluator : public StatementEvaluatorInterface<_Real> {

      using Real = DECLARE_DEFAULT(_Real, double);

    public:

      /*******************************************************************************
       * Section: Start of interface definition
       *
       */
      using Handle = PrimalTapeStatementFunctions const*;

      template<typename Tape, typename ... Args>
      static Real callForward(Handle const& h, Args&& ... args) {
        return ((FunctionForward<Tape>)h->forward)(std::forward<Args>(args)...);
      }

      template<typename Tape, typename ... Args>
      static Real callPrimal(Handle const& h, Args&& ... args) {
        return ((FunctionPrimal<Tape>)h->primal)(std::forward<Args>(args)...);
      }

      template<typename Tape, typename ... Args>
      static void callReverse(Handle const& h, Args&& ... args) {
        ((FunctionReverse<Tape>)h->reverse)(std::forward<Args>(args)...);
      }

      template<typename Tape, typename Generator, typename Expr>
      static Handle createHandle() {
        return &DirectStatementEvaluatorStaticStore<Generator, Expr>::staticStore;
      }

    protected:

      template<typename Tape>
      using FunctionForward = decltype(&Tape::template statementEvaluateForward<ActiveType<Tape>>);

      template<typename Tape>
      using FunctionPrimal = decltype(&Tape::template statementEvaluatePrimal<ActiveType<Tape>>);

      template<typename Tape>
      using FunctionReverse = decltype(&Tape::template statementEvaluateReverse<ActiveType<Tape>>);


  };
}
