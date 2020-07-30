#pragma once

#include "../../aux/macros.hpp"
#include "../../config.h"

/** \copydoc codi::Namespace */
namespace codi {

  template<typename _Real, typename _Gradient, typename _Identifier>
  struct ManualStatementPushTapeInterface {
    public:

      using Real = CODI_DECLARE_DEFAULT(_Real, double);
      using Gradient = CODI_DECLARE_DEFAULT(_Gradient, double);
      using Identifier = CODI_DECLARE_DEFAULT(_Identifier, int);

      /*******************************************************************************
       * Section: Start of interface definition
       *
       */

      void pushJacobiManual(Real const& jacobi, Real const& value, Identifier const& index);
      void storeManual(Real const& lhsValue, Identifier& lhsIndex, Config::ArgumentSize const& size);

  };
}
