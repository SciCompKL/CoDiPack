#pragma once

#include "../../aux/macros.hpp"
#include "../../config.h"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Add derivative information for custom operations on the tape.
   *
   * See \ref TapeInterfaces for a general overview of the tape interface design in CoDiPack.
   *
   * The function in this interface can be used to provide derivative information to CoDiPack for functions, that are
   * to small to be evaluated via an external function.
   *
   * The \ref sec_forwardAD "forward" and \ref sec_reverseAD "reverse" AD equations are the base for this interface. The user has
   * to provide the Jacobian \f$ \frac{\d phi}{\d u} \f$ for all arguments \f$u\f$ and compute the value for \f$w\f$.
   *
   * Before the call to storeManual the user has to update the value of the output, that is \f$w\f$ in the above
   * equations. This is usually done with `output.value() = w`. Afterwards storeManual() has to be called. The `size`
   * argument is the number of arguments \f$u\f$ from the equations above.
   *
   * Afterwards the user has to call pushJacobiManual() for each argument \f$u\f$.
   *
   * The user has to ensure, that the computations for the Jacobians are evaluated such that the CoDiPack tape does not
   * recognize them.
   *
   * An example manual push is(documentation/examples/manualStatementPushTapeInterface.cpp):
   * \snippet examples/manualStatementPushTapeInterface.cpp Manual statement push
   *
   * @tparam _Real        The computation type of a tape usually defined by ActiveType::Real.
   * @tparam _Gradient    The gradient type of a tape usually defined by ActiveType::Gradient.
   * @tparam _Identifier  The adjoint/tangent identification of a tape usually defined by ActiveType::Identifier.
   */
  template<typename _Real, typename _Gradient, typename _Identifier>
  struct ManualStatementPushTapeInterface {
    public:

      using Real = CODI_DD(_Real, double); ///< See ManualStatementPushTapeInterface.
      using Gradient = CODI_DD(_Gradient, double); ///< See ManualStatementPushTapeInterface.
      using Identifier = CODI_DD(_Identifier, int); ///< See ManualStatementPushTapeInterface.

      /*******************************************************************************/
      /// @name Interface definition

       /// Push a Jacobian entry to the tape. storeManual() has to be called before. This method needs to be called
       /// as often as the size argument.
      void pushJacobiManual(Real const& jacobi, Real const& value, Identifier const& index);

      /// Initialize the storing of hand computed statement. The primal value needs to be already updated.
      void storeManual(Real const& lhsValue, Identifier& lhsIndex, Config::ArgumentSize const& size);

  };
}
