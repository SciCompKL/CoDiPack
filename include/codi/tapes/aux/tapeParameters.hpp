#pragma once

#include "../../config.h"

/** \copydoc codi::Namespace */
namespace codi {

  enum class TapeParameters {
    AdjointSize,
    PrimalSize,
    StatementSize,
    JacobianSize,
    PassiveValuesSize,
    ConstantValuesSize,
    RhsIdentifiersSize,
    ExternalFunctionsSize
  };
}
