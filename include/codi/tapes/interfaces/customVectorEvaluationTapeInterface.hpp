#pragma once

#include "../../aux/macros.hpp"
#include "../../config.h"
#include "../data/position.hpp"
#include "forwardEvaluationTapeInterface.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  template<typename _Position>
  struct CustomVectorEvaluationTapeInterface : public virtual ForwardEvaluationTapeInterface<_Position> {
    public:

      using Position = CODI_DECLARE_DEFAULT(_Position, EmptyPosition);

      /*******************************************************************************
       * Section: Start of interface definition
       *
       */

      template<typename Adjoint>
      void evaluate(Position const& start, Position const& end, Adjoint* data);
      template<typename Adjoint>
      void evaluateForward(Position const& start, Position const& end, Adjoint* data);
  };
}
