#pragma once

#include "../../aux/macros.h"
#include "../../config.h"
#include "../data/position.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  template<typename _Position>
  struct PositionalEvaluationTapeInterface {
    public:

      using Position = DECLARE_DEFAULT(_Position, EmptyPosition);

      /*******************************************************************************
       * Section: Start of interface definition
       *
       */

      void evaluate(Position const& start, Position const& end);
      void clearAdjoints(Position const& start, Position const& end);

      Position getPosition() const;
      Position getZeroPosition() const;

      void resetTo(Position const& pos);
  };
}
