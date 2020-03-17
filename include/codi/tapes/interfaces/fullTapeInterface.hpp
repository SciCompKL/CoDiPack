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
      public CustomVectorEvaluationTapeInterface<_Position>,
      public DataManagementTapeInterface,
      public ExternalFunctionTapeInterface<_Real, _Gradient, _Identifier>,
      public ForwardEvaluationTapeInterface<_Position>,
      public GradientAccessTapeInterface<_Gradient, _Identifier>,
      public IdentifierInformationTapeInterface<_Real, _Gradient, _Identifier>,
      public InternalExpressionTapeInterface<_Identifier>,
      public ManualStatementPushTapeInterface<_Real, _Gradient, _Identifier>,
      public PositionalEvaluationTapeInterface<_Position>,
      public PreaccumulationEvaluationTapeInterface<_Real, _Gradient, _Identifier, _Position>,
      public PrimalEvaluationTapeInterface<_Real, _Identifier, _Position>,
      public ReverseTapeInterface<_Real, _Gradient, _Identifier>
  {
    public:

      using Real = DECLARE_DEFAULT(_Real, double);
      using Gradient = DECLARE_DEFAULT(_Gradient, double);
      using Identifier = DECLARE_DEFAULT(_Identifier, int);
      using Position = DECLARE_DEFAULT(_Position, EmptyPosition);

  };
}
