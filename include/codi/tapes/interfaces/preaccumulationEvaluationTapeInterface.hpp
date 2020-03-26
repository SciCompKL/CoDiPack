#pragma once
#pragma once

#include "../../aux/macros.h"
#include "../../config.h"
#include "../data/position.hpp"
#include "positionalEvaluationTapeInterface.hpp"
#include "forwardEvaluationTapeInterface.hpp"
#include "manualStatementPushTapeInterface.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  template<typename _Real, typename _Gradient, typename _Identifier, typename _Position>
  struct PreaccumulationEvaluationTapeInterface :
      public virtual PositionalEvaluationTapeInterface<_Position>,
      public virtual ForwardEvaluationTapeInterface<_Position>,
      public virtual ManualStatementPushTapeInterface<_Real, _Gradient, _Identifier>
  {
    public:

      using Real = DECLARE_DEFAULT(_Real, double);
      using Gradient = DECLARE_DEFAULT(_Gradient, double);
      using Identifier = DECLARE_DEFAULT(_Identifier, int);
      using Position = DECLARE_DEFAULT(_Position, EmptyPosition);

      /*******************************************************************************
       * Section: Start of interface definition
       *
       */

      void evaluateKeepState(Position const& start, Position const& end);
      void evaluateForwardKeepState(Position const& start, Position const& end);

  };
}
