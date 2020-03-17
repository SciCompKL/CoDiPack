#pragma once

#include "../../aux/macros.h"
#include "../../config.h"
#include "../../expressions/lhsExpressionInterface.hpp"
#include "../aux/externalFunction.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  template<typename _Real, typename _Gradient, typename _Identifier>
  struct ExternalFunctionTapeInterface {
    public:

      using Real = DECLARE_DEFAULT(_Real, double);
      using Gradient = DECLARE_DEFAULT(_Gradient, double);
      using Identifier = DECLARE_DEFAULT(_Identifier, int);

      /*******************************************************************************
       * Section: Start of interface definition
       *
       */

      template<typename Lhs>
      Real registerExternalFunctionOutput(LhsExpressionInterface<Real, Gradient, CustomVectorEvaluationTapeInterface, Lhs>& value);

      void pushExternalFunction(ExternalFunction const& extFunc);
      // TODO: add typed function push
  };
}
