#pragma once

#include "../../config.h"

/** \copydoc codi::Namespace */
namespace codi {

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
