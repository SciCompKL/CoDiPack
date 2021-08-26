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
   * The definitions in the ReverseTapeInterface provide only methods that operate on the full tape. With the methods
   * in this interface, all these operations can be performed on parts of the tape.
   *
   * Here is an example for a positional tape evaluation (\ref Example_03_Positional_tape_evaluations):
   * \snippet examples/Example_03_Positional_tape_evaluations.cpp Positional evaluation
   *
   * @tparam T_Position  Global tape position, usually chosen as Tape::Position.
   */
  template<typename T_Position>
  struct PositionalEvaluationTapeInterface {
    public:

      using Position = CODI_DD(T_Position, EmptyPosition);  ///< See PositionalEvaluationTapeInterface.

      /*******************************************************************************/
      /// @name Interface definition

      /// Perform a reverse evaluation for a part of the tape. It hast to hold start >= end.
      void evaluate(Position const& start, Position const& end);

      /// Clear all adjoints that would be set in a tape evaluation from start to end. It has to hold start >= end.
      void clearAdjoints(Position const& start, Position const& end);

      Position getPosition() const;      ///< Current position of the tape.
      Position getZeroPosition() const;  ///< Initial position of the tape.

      void resetTo(Position const& pos, bool resetAdjoints = true);  ///< Rest the tape to the provided position.
  };
}
