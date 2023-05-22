/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2023 Chair for Scientific Computing (SciComp), University of Kaiserslautern-Landau
 * Homepage: http://www.scicomp.uni-kl.de
 * Contact:  Prof. Nicolas R. Gauger (codi@scicomp.uni-kl.de)
 *
 * Lead developers: Max Sagebaum, Johannes Blühdorn (SciComp, University of Kaiserslautern-Landau)
 *
 * This file is part of CoDiPack (http://www.scicomp.uni-kl.de/software/codi).
 *
 * CoDiPack is free software: you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * CoDiPack is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty
 * of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 * See the GNU General Public License for more details.
 * You should have received a copy of the GNU
 * General Public License along with CoDiPack.
 * If not, see <http://www.gnu.org/licenses/>.
 *
 * For other licensing options please contact us.
 *
 * Authors:
 *  - SciComp, University of Kaiserslautern-Landau:
 *    - Max Sagebaum
 *    - Johannes Blühdorn
 *    - Former members:
 *      - Tim Albring
 */
#pragma once

#include "../../config.h"
#include "../../misc/macros.hpp"
#include "../data/position.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  template<typename T_Real, typename T_Gradient, typename T_Tape, typename T_Impl>
  struct LhsExpressionInterface;

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
   * @tparam T_Real        The computation type of a tape, usually chosen as ActiveType::Real.
   * @tparam T_Gradient    The gradient type of a tape, usually chosen as ActiveType::Gradient.
   * @tparam T_Identifier  The adjoint/tangent identification type of a tape, usually chosen as ActiveType::Identifier.
   */
  template<typename T_Real, typename T_Gradient, typename T_Identifier>
  struct IdentifierInformationTapeInterface {
    public:

      using Real = CODI_DD(T_Real, double);           ///< See IdentifierInformationTapeInterface.
      using Gradient = CODI_DD(T_Gradient, double);   ///< See IdentifierInformationTapeInterface.
      using Identifier = CODI_DD(T_Identifier, int);  ///< See IdentifierInformationTapeInterface.

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
      ///
      /// @tparam Lhs  Class that implements the LhsExpressionInterface. See also LhsExpressionInterface.
      /// @tparam Tape  Tape implementation used in the LhsExpressionInterface. See also LhsExpressionInterface.
      template<typename Lhs, typename Tape>
      void deactivateValue(LhsExpressionInterface<Real, Gradient, Tape, Lhs>& value);
  };
}
