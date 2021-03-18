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
    AdjointSize,         ///< [A: RW] Current adjoint vector size, not the maximum possible size. See LargestIdentifier.
    ConstantValuesSize,  ///< [A: RW] Allocated size of the constant value vector in primal value tapes.
    ExternalFunctionsSize,  ///< [A: RW] Allocated size of the external function vector.
    JacobianSize,           ///< [A: RW] Allocated size of the argument Jacobian vector in Jacobian tapes.
    LargestIdentifier,      ///< [A: R] Largest identifier distributed by the index manger.
    PassiveValuesSize,      ///< [A: RW] Allocated size of the passive value vector in primal value tapes.
    PrimalSize,             ///< [A: RW] Primal vector size in primal value tapes.
    RhsIdentifiersSize,     ///< [A: RW] Allocated size of the right hand side identifiers vector in primal value tapes.
    StatementSize           ///< [A: RW] Allocated size of the statement vector in all tapes.
  };
}
