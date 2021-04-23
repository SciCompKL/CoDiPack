#pragma once

#include "../../aux/macros.hpp"
#include "../../config.h"
#include "../data/position.hpp"
#include "positionalEvaluationTapeInterface.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Forward AD evaluation of a recorded tape.
   *
   * See \ref TapeInterfaces for a general overview of the tape interface design in CoDiPack.
   *
   * Here is an example for a forward tape evaluation (documentation/examples/forwardEvaluationTapeInterface.cpp):
   * \snippet examples/forwardEvaluationTapeInterface.cpp Forward tape evaluation
   *
   * @tparam _Position  Global tape position type, usually chosen as Tape::Position.
   */
  template<typename _Position>
  struct ForwardEvaluationTapeInterface : public virtual PositionalEvaluationTapeInterface<_Position> {
    public:

      using Position = CODI_DD(_Position, EmptyPosition);  ///< See ForwardEvaluationTapeInterface.

      /*******************************************************************************/
      /// @name Interface definition

      /// Perform a forward evaluation of a part of the tape. It hast to hold start <= end.
      void evaluateForward(Position const& start, Position const& end);

      /// Perform a forward evaluation of the full tape.
      void evaluateForward();
  };
}
