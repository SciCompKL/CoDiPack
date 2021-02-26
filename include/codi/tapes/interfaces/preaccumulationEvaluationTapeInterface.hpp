#pragma once
#pragma once

#include "../../aux/macros.hpp"
#include "../../config.h"
#include "../data/position.hpp"
#include "positionalEvaluationTapeInterface.hpp"
#include "forwardEvaluationTapeInterface.hpp"
#include "manualStatementPushTapeInterface.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Perform tape evaluations but ensure that the state does not change when it is finished.
   *
   * See \ref TapeInterfaces for a general overview of the tape interface design in CoDiPack.
   *
   * These interface functions are used for small tape evaluation where only a part of the tape is evaluated. Especially
   * for primal value tapes, it is essential that the primal value store in the background is in sync with the program
   * state. The normal evaluate methods in these tapes copy the primal value vector an perform all operations on the
   * copied vector. In the *KeepState() methods, they ensure the correctness of the primal value vector by performing
   * e.g. a primal evaluation after the reverse evaluation. This will yield better performance for small tape ranges in
   * the evaluation.
   *
   * @tparam _Real        The computation type of a tape usually defined by ActiveType::Real.
   * @tparam _Gradient    The gradient type of a tape usually defined by ActiveType::Gradient.
   * @tparam _Identifier  The adjoint/tangent identification of a tape usually defined by ActiveType::Identifier.
   * @tparam _Position  Global tape position usually defined by Tape::Position.
   */
  template<typename _Real, typename _Gradient, typename _Identifier, typename _Position>
  struct PreaccumulationEvaluationTapeInterface :
      public virtual PositionalEvaluationTapeInterface<_Position>,
      public virtual ForwardEvaluationTapeInterface<_Position>,
      public virtual ManualStatementPushTapeInterface<_Real, _Gradient, _Identifier>
  {
    public:

      using Real = CODI_DD(_Real, double);  ///< See PreaccumulationEvaluationTapeInterface.
      using Gradient = CODI_DD(_Gradient, double);  ///< See PreaccumulationEvaluationTapeInterface.
      using Identifier = CODI_DD(_Identifier, int);  ///< See PreaccumulationEvaluationTapeInterface.
      using Position = CODI_DD(_Position, EmptyPosition);  ///< See PreaccumulationEvaluationTapeInterface.

      /*******************************************************************************/
      /// @name Interface definition

      /// Perform a tape evaluation but restore the state afterwards such that it is the same as when the evaluation started.
      void evaluateKeepState(Position const& start, Position const& end);

      /// Perform a tape evaluation but restore the state afterwards such that it is the same as when the evaluation started.
      void evaluateForwardKeepState(Position const& start, Position const& end);

  };
}
