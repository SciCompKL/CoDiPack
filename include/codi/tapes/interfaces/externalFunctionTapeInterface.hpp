#pragma once

#include "../../aux/macros.hpp"
#include "../../config.h"
#include "../../expressions/lhsExpressionInterface.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  template<typename _Tape>
  struct ExternalFunction;

  /**
   * @brief Add user defined functions to the tape evaluation.
   *
   * See \ref TapeInterfaces for a general overview of the tape interface design in CoDiPack.
   *
   * External functions allow the user to evaluate custom operations during a tape evaluation. Each external function
   * has pointers for the reverse, forward and primal evaluation of a tape. A function pointer may be null if the
   * corresponding mode is not called on the tape. Otherwise, if the corresponding component is not defined in the
   * external function, then a CODI_EXCEPTION is thrown.
   *
   * What kind of operations are evaluated in the external function is up to the user. They are usually used to define
   * derivative computations for libraries that cannot be differentiated with operator overloading.
   *
   * Variables that are outputs of external functions have to be registered with registerExternalFunctionOutput. This
   * will ensure that the variable is considered as active in CoDiPack. For primal value tapes, the return value of this
   * function provides the old value stored under the identifier that the variable has received. This old value has to
   * be restored with a call to adjointInterface.setPrimal() during the evaluation of the external function in reverse
   * mode.
   *
   * Here is an example (documentation/examples/externalFunctionTapeInterface.cpp):
   * \snippet examples/externalFunctionTapeInterface.cpp External function
   *
   * @tparam _Real        The computation type of a tape, usually chosen as ActiveType::Real.
   * @tparam _Gradient    The gradient type of a tape, usually chosen as ActiveType::Gradient.
   * @tparam _Identifier  The adjoint/tangent identification type of a tape, usually chosen as ActiveType::Identifier.
   */
  template<typename _Real, typename _Gradient, typename _Identifier>
  struct ExternalFunctionTapeInterface {
    public:

      using Real = CODI_DD(_Real, double);           ///< See ExternalFunctionTapeInterface.
      using Gradient = CODI_DD(_Gradient, double);   ///< See ExternalFunctionTapeInterface.
      using Identifier = CODI_DD(_Identifier, int);  ///< See ExternalFunctionTapeInterface.

      /*******************************************************************************/
      /// @name Interface definition

      /// Register an external function output on the tape.
      /// @return For primal value tapes, the return value has to be stored by the external function. The value has to
      ///         be restored with a call to adjointInterface.setPrimal() during the evaluation of the external function
      ///         in reverse mode. For this purpose, the primal value is identified by the index which the variable
      ///         received when it was registered with registerExternalFunctionOutput.
      template<typename Lhs>
      Real registerExternalFunctionOutput(
          LhsExpressionInterface<Real, Gradient, ExternalFunctionTapeInterface, Lhs>& value);

      /// Push an external function to the tape.
      void pushExternalFunction(ExternalFunction<ExternalFunctionTapeInterface> const& extFunc);
  };
}
