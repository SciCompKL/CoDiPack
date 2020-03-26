#pragma once

#include "../../aux/macros.h"
#include "../../config.h"

#include "customVectorEvaluationTapeInterface.hpp"
#include "dataManagementTapeInterface.hpp"
#include "externalFunctionTapeInterface.hpp"
#include "forwardEvaluationTapeInterface.hpp"
#include "gradientAccessTapeInterface.hpp"
#include "identifierInformationTapeInterface.hpp"
#include "internalExpressionTapeInterface.hpp"
#include "manualStatementPushTapeInterface.hpp"
#include "positionalEvaluationTapeInterface.hpp"
#include "preaccumulationEvaluationTapeInterface.hpp"
#include "primalEvaluationTapeInterface.hpp"
#include "reverseTapeInterface.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  template<typename _Real, typename _Gradient, typename _Identifier, typename _Position>
  struct FullTapeInterface :
      public virtual CustomVectorEvaluationTapeInterface<_Position>,
      public virtual DataManagementTapeInterface,
      public virtual ExternalFunctionTapeInterface<_Real, _Gradient, _Identifier>,
      public virtual ForwardEvaluationTapeInterface<_Position>,
      public virtual GradientAccessTapeInterface<_Gradient, _Identifier>,
      public virtual IdentifierInformationTapeInterface<_Real, _Gradient, _Identifier>,
      public virtual InternalExpressionTapeInterface<_Identifier>,
      public virtual ManualStatementPushTapeInterface<_Real, _Gradient, _Identifier>,
      public virtual PositionalEvaluationTapeInterface<_Position>,
      public virtual PreaccumulationEvaluationTapeInterface<_Real, _Gradient, _Identifier, _Position>,
      public virtual PrimalEvaluationTapeInterface<_Real, _Identifier, _Position>,
      public virtual ReverseTapeInterface<_Real, _Gradient, _Identifier>
  {
    public:

      using Real = DECLARE_DEFAULT(_Real, double);
      using Gradient = DECLARE_DEFAULT(_Gradient, double);
      using Identifier = DECLARE_DEFAULT(_Identifier, int);
      using Position = DECLARE_DEFAULT(_Position, EmptyPosition);

  };
}
