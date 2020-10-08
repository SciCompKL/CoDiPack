#pragma once

#include "../../aux/macros.hpp"
#include "../../config.h"
#include "../data/position.hpp"
#include "positionalEvaluationTapeInterface.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Forward AD evaluation of recorded tape.
   *
   * See \ref TapeInterfaces for a general overview of the tape interface design in CoDiPack.
   *
   * An example for a forward tape evaluation (documentation/examples/forwardEvaluationTapeInterface.cpp):
   * \snippet examples/forwardEvaluationTapeInterface.cpp Forward tape evaluation
   *
   * @tparam _Position  Global tape position usually defined by Tape::Position.
   */
  template<typename _Position>
  struct ForwardEvaluationTapeInterface : public virtual PositionalEvaluationTapeInterface<_Position> {
    public:

      using Position = CODI_DECLARE_DEFAULT(_Position, EmptyPosition); ///< See ForwardEvaluationTapeInterface.

      /*******************************************************************************/
      /// @name Interface definition

      /// Perform a forward evaluation for a part of the tape. It hast to hold start <= end.
      void evaluateForward(Position const& start, Position const& end);

      /// Perform a forward evaluation for the full.
      void evaluateForward();
  };
}
