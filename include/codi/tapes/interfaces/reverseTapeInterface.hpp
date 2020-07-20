#pragma once

#include <iostream>

#include "../../aux/macros.h"
#include "../../config.h"
#include "../../expressions/lhsExpressionInterface.hpp"
#include "../../expressions/logic/traversalLogic.hpp"
#include "../../traits/expressionTraits.hpp"
#include "../aux/tapeValues.hpp"
#include "gradientAccessTapeInterface.hpp"
#include "internalStatementRecordingInterface.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  template<typename _Real, typename _Gradient, typename _Identifier>
  struct ReverseTapeInterface : public virtual InternalStatementRecordingInterface<_Identifier>,
                                public virtual GradientAccessTapeInterface<_Gradient, _Gradient> {


      using Real = DECLARE_DEFAULT(_Real, double);
      using Gradient = DECLARE_DEFAULT(_Gradient, double);
      using Identifier = DECLARE_DEFAULT(_Identifier, int);

      /*******************************************************************************
       * Section: Start of interface definition
       *
       */

      void evaluate();

      template<typename Lhs> void registerInput(LhsExpressionInterface<Real, Gradient, ReverseTapeInterface, Lhs>& value);
      template<typename Lhs> void registerOutput(LhsExpressionInterface<Real, Gradient, ReverseTapeInterface, Lhs>& value);

      void setActive();
      void setPassive();
      bool isActive() const;

      void clearAdjoints();
      void reset(bool resetAdjoints = true);

      template<typename Stream = std::ostream> void printStatistics(Stream& out = std::cout) const;
      template<typename Stream = std::ostream> void printTableHeader(Stream& out = std::cout) const;
      template<typename Stream = std::ostream> void printTableRow(Stream& out = std::cout) const;
      TapeValues getTapeValues() const;
  };
}

