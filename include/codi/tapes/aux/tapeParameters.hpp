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
