#pragma once

#include "../../aux/macros.hpp"
#include "../../config.h"
#include "../../expressions/lhsExpressionInterface.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  template<typename _Tape>
  struct ExternalFunction;

  /**
   * @brief Add user defined function to the tape evaluation.
   *
   * See \ref TapeInterfaces for a general overview of the tape interface design in CoDiPack.
   *
   * External functions allow the user to evaluate custom operations during a tape evaluation. Each external function
   * has pointers for the reverse, forward and primal evaluation of an tape. The function pointer only needs to be not
   * null if the corresponding mode is called on the tape. If the proper function is not defined in an external function,
   * then a CODI_EXCEPTION is thrown.
   *
   * What kind of operations are evaluated in the external function are up to the user. Usually they are used to define
   * derivative computations for libraries, that can not be differentiated with operator overloading.
   *
   * Variables that are outputs of external functions, need to be registered with registerExternalFunctionOutput. This
   * will ensure that the variable is considered as active in CoDiPack. For primal value tapes the return value of this
   * function provides the old value stored under this identifier, that the value has received. This old value needs to
   * be restored with a call to adjointInterface.setPrimal() during the evaluation of the external function in a reverse
   * mode.
   *
   * An example is (documentation/examples/externalFunctionTapeInterface.cpp):
   * \snippet examples/externalFunctionTapeInterface.cpp External function
   *
   * @tparam _Real        The computation type of a tape usually defined by ActiveType::Real.
   * @tparam _Gradient    The gradient type of a tape usually defined by ActiveType::Gradient.
   * @tparam _Identifier  The adjoint/tangent identification of a tape usually defined by ActiveType::Identifier.
   */
  template<typename _Real, typename _Gradient, typename _Identifier>
  struct ExternalFunctionTapeInterface {
    public:

      using Real = CODI_DD(_Real, double); ///< See ExternalFunctionTapeInterface
      using Gradient = CODI_DD(_Gradient, double); ///< See ExternalFunctionTapeInterface
      using Identifier = CODI_DD(_Identifier, int); ///< See ExternalFunctionTapeInterface

      /*******************************************************************************/
      /// @name Interface definition


      /// Register an external function output on the tape.
      /// @return For primal value tapes the return value needs to be stored by the external function. The value needs
      ///         to be restored with a call to adjointInterface.setPrimal(). The index is the index value has gained by
      ///         registering it with this method.
      template<typename Lhs>
      Real registerExternalFunctionOutput(LhsExpressionInterface<Real, Gradient, ExternalFunctionTapeInterface, Lhs>& value);

       /// Push an external function to the tape.
      void pushExternalFunction(ExternalFunction<ExternalFunctionTapeInterface> const& extFunc);
  };
}
