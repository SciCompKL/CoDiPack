#pragma once

#include "../../config.h"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Possible configuration options for a tape.
   *
   * See DataManagementTapeInterface for details.
   */
  enum class TapeParameters {
    AdjointSize,
    ConstantValuesSize,
    ExternalFunctionsSize,
    JacobianSize,
    LargestIdentifier,
    PassiveValuesSize,
    PrimalSize,
    RhsIdentifiersSize,
    StatementSize
  };
}
