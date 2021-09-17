/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2021 Chair for Scientific Computing (SciComp), TU Kaiserslautern
 * Homepage: http://www.scicomp.uni-kl.de
 * Contact:  Prof. Nicolas R. Gauger (codi@scicomp.uni-kl.de)
 *
 * Lead developers: Max Sagebaum, Johannes Blühdorn (SciComp, TU Kaiserslautern)
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
 *  - SciComp, TU Kaiserslautern:
 *    - Max Sagebaum
 *    - Johannes Blühdorn
 *    - Former members:
 *      - Tim Albring
 */
#pragma once

#include "../../config.h"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Configuration options for a tape.
   *
   * See DataManagementTapeInterface for details.
   *
   * Access is defined by [A: "access"]. Options are:
   *  - R read only access (getParameter)
   *  - W write only access (setParameter)
   *  - RW read and write access (getParameter and setParameters)
   */
  enum class TapeParameters
  {
    AdjointSize,         ///< [A: RW] Current number of adjoint vector entries, not the maximum possible size.
                         ///<         See LargestIdentifier.
    ConstantValuesSize,  ///< [A: RW] Allocated number of entries in  the constant value vector in primal value tapes.
    ExternalFunctionsSize,  ///< [A: RW] Allocated number of entries in the external function vector.
    JacobianSize,           ///< [A: RW] Allocated number of entries in the argument Jacobian vector in Jacobian tapes.
    LargestIdentifier,      ///< [A: R] Largest identifier distributed by the index manger.
    PassiveValuesSize,      ///< [A: RW] Allocated number of entries in the passive value vector in primal value tapes.
    PrimalSize,             ///< [A: RW] Number of primal vector entries in primal value tapes.
    RhsIdentifiersSize,     ///< [A: RW] Allocated number of entries in the right hand side identifiers vector in primal
                            ///<         value tapes.
    StatementSize           ///< [A: RW] Allocated number of entries in the statement vector in all tapes.
  };
}
