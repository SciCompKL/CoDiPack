#pragma once

#include "../../aux/macros.h"
#include "../../config.h"

/** \copydoc codi::Namespace */
namespace codi {

  template<typename _Real, typename _Gradient, typename _Identifier>
  struct ManualStatementPushTapeInterface {
    public:

      using Real = DECLARE_DEFAULT(_Real, double);
      using Gradient = DECLARE_DEFAULT(_Gradient, double);
      using Identifier = DECLARE_DEFAULT(_Identifier, int);

      /*******************************************************************************
       * Section: Start of interface definition
       *
       */

      void pushJacobiManual(Real const& jacobi, Real const& value, Gradient const& index);
      void storeManual(Real const& lhsValue, Gradient& lhsIndex, Config::ArgumentSize const& size);

  };
}
