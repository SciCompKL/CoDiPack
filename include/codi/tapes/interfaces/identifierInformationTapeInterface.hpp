#pragma once

#include "../../aux/macros.hpp"
#include "../../config.h"
#include "../../expressions/lhsExpressionInterface.hpp"
#include "../data/position.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief General information about the identifiers and checks if variables are active.
   *
   * See \ref TapeInterfaces for a general overview of the tape interface design in CoDiPack.
   *
   * With this interface, the user can check if a variable in the program is active or not. For an explanation of what
   * is an active variable for CoDiPack, please see \ref ActivityAnalysis.
   *
   * Here is an example for deactivating an identifier (documentation/examples/identifierInformationTapeInterface.cpp):
   * \snippet examples/identifierInformationTapeInterface.cpp Identifier Activity
   *
   * @tparam _Real        The computation type of a tape, usually chosen as ActiveType::Real.
   * @tparam _Gradient    The gradient type of a tape, usually chosen as ActiveType::Gradient.
   * @tparam _Identifier  The adjoint/tangent identification type of a tape, usually chosen as ActiveType::Identifier.
   */
  template<typename _Real, typename _Gradient, typename _Identifier>
  struct IdentifierInformationTapeInterface {
    public:

      using Real = CODI_DD(_Real, double);           ///< See IdentifierInformationTapeInterface.
      using Gradient = CODI_DD(_Gradient, double);   ///< See IdentifierInformationTapeInterface.
      using Identifier = CODI_DD(_Identifier, int);  ///< See IdentifierInformationTapeInterface.

      /*******************************************************************************/
      /// @name Interface definition

      /// True if the tape uses an index handler that provides identifiers in a monotonically increasing way (see
      /// LinearIndexManager).
      static bool constexpr LinearIndexHandling = CODI_UNDEFINED_VALUE;

      Identifier getPassiveIndex() const;                      ///< Identifier for passive values. Usually 0.
      Identifier getInvalidIndex() const;                      ///< Invalid identifier.
      bool isIdentifierActive(Identifier const& index) const;  ///< True if the identifier is considered active by the
                                                               ///< tape.

      /// Modify the value such that it is no longer active.
      template<typename Lhs>
      void deactivateValue(LhsExpressionInterface<Real, Gradient, IdentifierInformationTapeInterface, Lhs>& value);
  };
}
