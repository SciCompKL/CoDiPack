#pragma once

#include "../../aux/macros.hpp"
#include "../../config.h"
#include "../data/position.hpp"
#include "positionalEvaluationTapeInterface.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Perform a primal reevaluation of the tape.
   *
   * See \ref TapeInterfaces for a general overview of the tape interface design in CoDiPack.
   *
   * If the tape manages primal values can be detected with static constant `HasPrimalValues`.
   *
   * In a primal value tape the correctness of the primal values is very important. The tapes should be programmed such
   * that the primal values stored in the tape are always up to date with the state of the program. Only through user
   * interaction this sync in states can be broken, but then the user should know what he is doing.
   *
   * The primal evaluation is used to reevaluate the internal primal of the tape at a different position. Since CoDiPack
   * does not treat if branches or other conditional elements in the code, the user has to ensure that this reevaluation
   * still provides the correct numerical computation.
   *
   * An example for a primal reevaluation (documentation/examples/primalEvaluationTapeInterface.cpp):
   * \snippet examples/primalEvaluationTapeInterface.cpp Primal evaluation
   *
   * @tparam _Real        The computation type of a tape usually defined by ActiveType::Real.
   * @tparam _Identifier  The adjoint/tangent identification of a tape usually defined by ActiveType::Identifier.
   * @tparam _Position  Global tape position usually defined by Tape::Position.
   */
  template<typename _Real, typename _Identifier, typename _Position>
  struct PrimalEvaluationTapeInterface : public virtual PositionalEvaluationTapeInterface<_Position> {
    public:

      using Real = CODI_DECLARE_DEFAULT(_Real, double); ///< See PrimalEvaluationTapeInterface
      using Identifier = CODI_DECLARE_DEFAULT(_Identifier, int); ///< See PrimalEvaluationTapeInterface
      using Position = CODI_DECLARE_DEFAULT(_Position, EmptyPosition); ///< See PrimalEvaluationTapeInterface

      /*******************************************************************************/
      /// @name Interface definition

      static bool constexpr HasPrimalValues = CODI_UNDEFINED_VALUE; ///< True if the tape has primal values
      static bool constexpr RequiresPrimalRestore = CODI_UNDEFINED_VALUE;  ///< True if the primal state changes during a reverse or forward evaluation

      /// Perform a partly (forward) reevaluation of the primals in the tape. It has to hold start <= end.
      void evaluatePrimal(Position const& start, Position const& end);

      /// Perform a full (forward) reevaluation of the primals in the tape
      void evaluatePrimal();

      void setPrimal(Identifier const& identifier, Real const& gradient); ///< Set primal value.
      Real const& getPrimal(Identifier const& identifier) const; ///< Get primal value.

      Real& primal(Identifier const& identifier);  ///< Writable reference to primal value.
      Real const& primal(Identifier const& identifier) const; ///< Read only reference to primal value.

      /// Revert the primals to the state indicates by pos.
      void revertPrimals(Position const& pos);
  };
}
