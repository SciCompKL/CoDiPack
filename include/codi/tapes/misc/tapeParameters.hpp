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
  enum class TapeParameters {
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

  /**
   * @brief Policies for management of the tape's interal adjoints.
   *
   * For the convenience of the user, tapes manage their internal adjoints automatically, which involves multiple
   * tasks. AdjointsManagement::Manual indicates that non of these tasks is performed - they are the responsibility of
   * the caller instead. Functions that take an AdjointsManagement parameter default to AdjointsManagement::Automatic
   * and document the individual effects of AdjointsManagement::Manual. An overview over all possible effects is given
   * below.
   *
   * <b>Bounds checking.</b> The function accesses the adjoints. In the automatic mode, it checks whether the
   * adjoints are sufficiently large. If they are not, they might be <b>resized</b> or the
   * function might work on or return <b>dummy values</b>. To optimize the memory usage and/or reduce the number of
   * reallocations, AdjointsManagement::Manual can be used to skip bounds checking and resizing. It is the
   * responsibility of the caller to ensure sufficient adjoints size, for example by calls to
   * DataManagementTapeInterface::resizeAdjointVector.
   *
   * <b>Declaration of adjoints usage (locking).</b> If a tape implements it adjoints against
   * InternalAdjointsInterface, it keeps track of whether the adjoint vector is in use, which is for example the case
   * during tape evaluations. This is to ensure mutual exclusion with reallocations, this is particularly important in
   * shared-memory parallel taping, see also ThreadSafeGlobalAdjoints. Declaration of usage involves setting a lock,
   * which can become a bottleneck if it is done frequently. To optimize the performance, multiple operations can be
   * grouped into a single usage declaration, by surrounding them by manual
   * DataManagementTapeInterface::beginUseAdjoints and DataManagementTapeInterface::endUseAdjoints calls and invoking
   * them with AdjointsManagement::Manual. Note that any method that results in adjoint vector resizing must be called
   * outside usage declarations, otherwise there would be a deadlock.
   */
  enum class AdjointsManagement {
    Manual,    ///< Do not perform any bounds checking, locking, or resizing.
    Automatic  ///< Manage internal adjoints automatically, including locking, bounds checking, and resizing.
  };
}
