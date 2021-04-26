#pragma once

#include "../../aux/macros.hpp"
#include "../../config.h"
#include "../data/position.hpp"
#include "forwardEvaluationTapeInterface.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Allows user defined vectors for the forward and adjoint evaluation.
   *
   * See \ref TapeInterfaces for a general overview of the tape interface design in CoDiPack.
   *
   * The two additional evaluate methods allow for the evaluation of the tape with a custom adjoint vector. The type
   * of the vector must support the following operators:
   *  - operator =
   *  - operator *(Tape::Real, Adjoint) (Scalar multiplication from the left)
   *  - operator +=
   * It must also specialize #codi::GradientTraits::TraitsImplementation
   *
   * Here is an example for an evaluation with a custom adjoint vector
   * (documentation/examples/customAdjointVectorEvaluationTapeInterface.cpp):
   * \snippet examples/customAdjointVectorEvaluationTapeInterface.cpp Custom vector
   *
   * @tparam _Position  Global tape position, usually chosen as Tape::Position.
   */
  template<typename _Position>
  struct CustomAdjointVectorEvaluationTapeInterface : public virtual ForwardEvaluationTapeInterface<_Position> {
    public:

      using Position = CODI_DD(_Position, EmptyPosition);  ///< See CustomAdjointVectorEvaluationTapeInterface.

      /*******************************************************************************/
      /// @name Interface definition

      /**
       * \copydoc codi::PositionalEvaluationTapeInterface::evaluate
       *
       * @tparam Adjoint  See CustomAdjointVectorEvaluationTapeInterface documentation.
       */
      template<typename Adjoint>
      void evaluate(Position const& start, Position const& end, Adjoint* data);

      /**
       * \copydoc codi::ForwardEvaluationTapeInterface::evaluate(Position const&, Position const&)
       *
       * @tparam Adjoint  See CustomAdjointVectorEvaluationTapeInterface documentation
       */
      template<typename Adjoint>
      void evaluateForward(Position const& start, Position const& end, Adjoint* data);
  };
}
