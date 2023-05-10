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

  /**
   * @brief Add derivative information for custom operations to the tape.
   *
   * See \ref TapeInterfaces for a general overview of the tape interface design in CoDiPack.
   *
   * The functions in this interface can be used to provide derivative information to CoDiPack for functions that are
   * not known to CoDiPack but so small that an external function implementation is an overkill.
   *
   * The \ref sec_forwardAD "forward" and \ref sec_reverseAD "reverse" AD equations are the base for this interface. The
   * user has to provide the Jacobian \f$ \frac{\d \phi}{\d u} \f$ for all arguments \f$u\f$ and compute the value for
   * \f$w\f$.
   *
   * Before the call to storeManual the user has to update the value of the output, that is, \f$w\f$ in the above
   * equations. This is usually done with `output.value() = w`. Afterwards, storeManual() has to be called. The `size`
   * argument is the number of arguments \f$u\f$ from the equations above. This call ensures that `output` gets a
   * proper identifier and the dependency chain is not broken or wrong for this variable.
   *
   * Afterwards the user has to call pushJacobianManual() for each argument \f$u\f$.
   *
   * The user has to ensure that the computations of the Jacobians are evaluated such that the CoDiPack tape does not
   * accidentally record them.
   *
   * Here is an example for manual statement push (documentation/examples/manualStatementPushTapeInterface.cpp):
   * \snippet examples/manualStatementPushTapeInterface.cpp Manual statement push
   *
   * @tparam T_Real        The computation type of a tape, usually chosen as ActiveType::Real.
   * @tparam T_Gradient    The gradient type of a tape, usually chosen as ActiveType::Gradient.
   * @tparam T_Identifier  The adjoint/tangent identification type of a tape, usually chosen as ActiveType::Identifier.
   */
  template<typename T_Real, typename T_Gradient, typename T_Identifier>
  struct ManualStatementPushTapeInterface {
    public:

      using Real = CODI_DD(T_Real, double);           ///< See ManualStatementPushTapeInterface.
      using Gradient = CODI_DD(T_Gradient, double);   ///< See ManualStatementPushTapeInterface.
      using Identifier = CODI_DD(T_Identifier, int);  ///< See ManualStatementPushTapeInterface.

      /*******************************************************************************/
      /// @name Interface definition

      /// Push a Jacobian entry to the tape. storeManual() has to be called first and is passed the number of arguments.
      /// Afterwards, this method has to be called once for each argument.
      ///
      /// @param jacobian  Jacobian \f$ \frac{\d \phi}{\d u_i} \f$ of the argument \f$ u_i \f$.
      /// @param value   Value of the argument \f$ u_i \f$. Usually `u_i.value()`.
      /// @param index   Identifier of the argument \f$ u_i \f$. Usually `u_i.identifier()`.
      void pushJacobianManual(Real const& jacobian, Real const& value, Identifier const& index);

      /// Initialize the storing of a hand computed statement. The primal value has to be updated already.
      /// @param lhsValue   Value of the result \f$ w \f$. Usually `w.value()`.
      /// @param lhsIndex   Identifier of the result \f$ w \f$. Usually `w.identifier()`.
      /// @param size       Number of arguments of \f$ \phi \f$.
      void storeManual(Real const& lhsValue, Identifier& lhsIndex, Config::ArgumentSize const& size);
  };
}
