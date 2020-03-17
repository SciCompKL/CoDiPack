#pragma once

#include "../../aux/macros.h"
#include "../../config.h"
#include "../data/position.hpp"
#include "positionalEvaluationTapeInterface.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  template<typename _Position>
  struct ForwardEvaluationTapeInterface : public PositionalEvaluationTapeInterface<_Position> {
    public:

      using Position = DECLARE_DEFAULT(_Position, EmptyPosition);

      /*******************************************************************************
       * Section: Start of interface definition
       *
       */

      void evaluateForward(Position const& start, Position const& end);
      void evaluateForward();
  };
}
