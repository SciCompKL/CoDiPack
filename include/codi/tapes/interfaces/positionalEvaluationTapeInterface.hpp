#pragma once

#include "../../aux/macros.hpp"
#include "../../config.h"
#include "../data/position.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Reverse AD evaluation for parts of a recorded tape.
   *
   * See \ref TapeInterfaces for a general overview of the tape interface design in CoDiPack.
   *
   * The definitions in the ReverseTapeInterface provide only methods, that operate on the full tape. With the methods
   * in this interface all these operations can be performed on parts of the tape.
   *
   * An example for a positional tape evaluation (documentation/examples/positionalEvaluationTapeInterface.cpp):
   * \snippet examples/positionalEvaluationTapeInterface.cpp Forward tape evaluation
   *
   * @tparam _Position  Global tape position usually defined by Tape::Position.
   */
  template<typename _Position>
  struct PositionalEvaluationTapeInterface {
    public:

      using Position = CODI_DECLARE_DEFAULT(_Position, EmptyPosition);

      /*******************************************************************************/
      /// @name Interface definition

      /// Perform a reverse evaluation for a part of the tape. It hast to hold start >= end.
      void evaluate(Position const& start, Position const& end);

      /// Clear all adjoints that would be set in a tape evaluation from start to end. It has to hold start >= end.
      void clearAdjoints(Position const& start, Position const& end);

      Position getPosition() const; ///< Current position of the tape.
      Position getZeroPosition() const; ///< Initial position of the tape.

      void resetTo(Position const& pos, bool resetAdjoints = true); ///< Rest the tape to the provided position.
  };
}
