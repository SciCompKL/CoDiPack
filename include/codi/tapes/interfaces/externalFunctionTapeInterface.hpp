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

/** \copydoc codi::Namespace */
namespace codi {

  template<typename T_Tape>
  struct ExternalFunction;

  template<typename T_Real, typename T_Gradient, typename T_Tape, typename T_Impl>
  struct LhsExpressionInterface;

  /**
   * @brief Add user defined functions to the tape evaluation.
   *
   * See \ref TapeInterfaces for a general overview of the tape interface design in CoDiPack.
   *
   * External functions allow the user to evaluate custom operations during a tape evaluation. Each external function
   * has pointers for the reverse, forward and primal evaluation of a tape. A function pointer may be null if the
   * corresponding mode is not called on the tape. Otherwise, if the corresponding pointer is null and the
   * corresponding mode is called on the tape, then a CODI_EXCEPTION is thrown.
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
   * @tparam T_Real        The computation type of a tape, usually chosen as ActiveType::Real.
   * @tparam T_Gradient    The gradient type of a tape, usually chosen as ActiveType::Gradient.
   * @tparam T_Identifier  The adjoint/tangent identification type of a tape, usually chosen as ActiveType::Identifier.
   */
  template<typename T_Real, typename T_Gradient, typename T_Identifier>
  struct ExternalFunctionTapeInterface {
    public:

      using Real = CODI_DD(T_Real, double);           ///< See ExternalFunctionTapeInterface.
      using Gradient = CODI_DD(T_Gradient, double);   ///< See ExternalFunctionTapeInterface.
      using Identifier = CODI_DD(T_Identifier, int);  ///< See ExternalFunctionTapeInterface.

      /*******************************************************************************/
      /// @name Interface definition

      /// Register an external function output on the tape.
      /// @return For primal value tapes, the return value has to be stored by the external function. The value has to
      ///         be restored with a call to adjointInterface.setPrimal() during the evaluation of the external function
      ///         in reverse mode. For this purpose, the primal value is identified by the index which the variable
      ///         received when it was registered with registerExternalFunctionOutput.
      ///
      /// @tparam Lhs  Class that implements the LhsExpressionInterface. See also LhsExpressionInterface.
      /// @tparam Tape  Tape implementation used in the LhsExpressionInterface. See also LhsExpressionInterface.
      template<typename Lhs, typename Tape>
      Real registerExternalFunctionOutput(LhsExpressionInterface<Real, Gradient, Tape, Lhs>& value);

      /// Push an external function to the tape.
      ///
      /// The external function class can be created via the helper ExternalFunction::create.
      void pushExternalFunction(ExternalFunction<ExternalFunctionTapeInterface> const& extFunc);
  };
}
