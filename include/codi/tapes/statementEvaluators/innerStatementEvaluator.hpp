#pragma once

#include <algorithm>
#include <functional>
#include <type_traits>

#include "../../aux/macros.h"
#include "../../expressions/activeType.hpp"
#include "../../traits/expressionTraits.hpp"
#include "statementEvaluatorInterface.hpp"
#include "directStatementEvaluator.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  struct InnerPrimalTapeStatementData : public PrimalTapeStatementFunctions {
    public:

      using Base = PrimalTapeStatementFunctions;

      size_t maxActiveArguments;
      size_t maxConstantArguments;

      InnerPrimalTapeStatementData(
          size_t maxActiveArguments, size_t maxConstantArguments,
          typename Base::Handle forward, typename Base::Handle primal, typename Base::Handle reverse) :
        Base(forward, primal, reverse),
        maxActiveArguments(maxActiveArguments),
        maxConstantArguments(maxConstantArguments) {}
  };

  template<typename Tape, typename Expr>
  struct InnerStatementEvaluatorStaticStore {
    public:

      static InnerPrimalTapeStatementData const staticStore;
  };

  template<typename Tape, typename Expr>
  InnerPrimalTapeStatementData const InnerStatementEvaluatorStaticStore<Tape, Expr>::staticStore(
      NumberOfActiveTypeArguments<Expr>::value,
      NumberOfConstantTypeArguments<Expr>::value,
      (typename PrimalTapeStatementFunctions::Handle)Tape::template statementEvaluateForwardInner<Expr>,
      (typename PrimalTapeStatementFunctions::Handle)Tape::template statementEvaluatePrimalInner<Expr>,
      (typename PrimalTapeStatementFunctions::Handle)Tape::template statementEvaluateReverseInner<Expr>);

  template<typename _Real>
  struct InnerStatementEvaluator : public StatementEvaluatorInterface<_Real> {
    public:

      using Real = DECLARE_DEFAULT(_Real, double);

      /*******************************************************************************
       * Section: Start of interface definition
       *
       */
      using Handle = InnerPrimalTapeStatementData const*;

      template<typename Tape, typename ... Args>
      static Real callForward(Handle const& h, Args&& ... args) {
        return Tape::statementEvaluateForwardFull(
              (FunctionForward<Tape>)h->forward, h->maxActiveArguments, h->maxConstantArguments,
              std::forward<Args>(args)...);
      }

      template<typename Tape, typename ... Args>
      static Real callPrimal(Handle const& h, Args&& ... args) {
        return Tape::statementEvaluatePrimalFull(
              (FunctionPrimal<Tape>)h->primal, h->maxActiveArguments, h->maxConstantArguments,
              std::forward<Args>(args)...);
      }

      template<typename Tape, typename ... Args>
      static void callReverse(Handle const& h, Args&& ... args) {
        Tape::statementEvaluateReverseFull(
              (FunctionReverse<Tape>)h->reverse, h->maxActiveArguments, h->maxConstantArguments,
              std::forward<Args>(args)...);
      }

      template<typename Tape, typename Expr>
      static Handle createHandle() {
        return &InnerStatementEvaluatorStaticStore<Tape, Expr>::staticStore;
      }

    protected:

      template<typename Tape>
      using FunctionForward = decltype(&Tape::template statementEvaluateForwardInner<ActiveType<Tape>>);

      template<typename Tape>
      using FunctionPrimal = decltype(&Tape::template statementEvaluatePrimalInner<ActiveType<Tape>>);

      template<typename Tape>
      using FunctionReverse = decltype(&Tape::template statementEvaluateReverseInner<ActiveType<Tape>>);


  };
}
