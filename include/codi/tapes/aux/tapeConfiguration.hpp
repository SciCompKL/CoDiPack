#pragma once

#include "../../config.h"

/** \copydoc codi::Namespace */
namespace codi {

  enum class ConfigurationOption {
    AdjointSize,
    PrimalSize,
    StatementSize,
    JacobianSize,
    PassiveValuesSize,
    ConstantValuesSize,
    RhsIdentifiersSize,
  };
}
