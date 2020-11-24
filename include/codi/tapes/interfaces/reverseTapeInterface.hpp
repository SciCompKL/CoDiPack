#pragma once

#include <iostream>

#include "../../aux/macros.hpp"
#include "../../config.h"
#include "../../expressions/lhsExpressionInterface.hpp"
#include "../../expressions/logic/traversalLogic.hpp"
#include "../../traits/expressionTraits.hpp"
#include "../../traits/realTraits.hpp"
#include "../aux/tapeValues.hpp"
#include "gradientAccessTapeInterface.hpp"
#include "internalStatementRecordingInterface.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Minimum tape interface for a working reverse tape implementation.
   *
   * See \ref TapeInterfaces for a general overview of the tape interface design in CoDiPack.
   *
   * For an example on how to uses this interface to perform an AD reverse mode recording and evaluation of a program
   * pleas see tutorials (TODO: ref).
   *
   * Implementation hints:
   * A tape should only record information if it is active. That is everything between a call to
   * setActive() and setPassive(). A call to setActive() is not considered to reset the tape in CoDiPack. This is only
   * done with a call to reset(). This allows the user to skip unnecessary parts in the application by setting the tape
   * to passive for these regions.
   *
   * A example use of a tape is(documentation/examples/reverseModeAD.cpp):
   * \snippet examples/reverseModeAD.cpp Reverse mode AD
   *
   * @tparam _Real        The computation type of a tape usually defined by ActiveType::Real.
   * @tparam _Gradient    The gradient type of a tape usually defined by ActiveType::Gradient.
   * @tparam _Identifier  The adjoint/tangent identification of a tape usually defined by ActiveType::Identifier.
   */
  template<typename _Real, typename _Gradient, typename _Identifier>
  struct ReverseTapeInterface : public virtual InternalStatementRecordingInterface<_Identifier>,
                                public virtual GradientAccessTapeInterface<_Gradient, _Gradient> {


      using Real = CODI_DECLARE_DEFAULT(_Real, double); ///< See ReverseTapeInterface
      using Gradient = CODI_DECLARE_DEFAULT(_Gradient, double); ///< See ReverseTapeInterface
      using Identifier = CODI_DECLARE_DEFAULT(_Identifier, int); ///< See ReverseTapeInterface

      using PassiveReal = PassiveRealType<Real>;

      /*******************************************************************************/
      /// @name Recording

      /// Mark a value as input (independent) and make it active
      template<typename Lhs> void registerInput(LhsExpressionInterface<Real, Gradient, ReverseTapeInterface, Lhs>& value);
      /// Mark a value as output (dependent)
      template<typename Lhs> void registerOutput(LhsExpressionInterface<Real, Gradient, ReverseTapeInterface, Lhs>& value);

      void setActive();      ///< Start/continue recording of statements
      void setPassive();     ///< Stop recording of statements
      bool isActive() const; ///< Check if the tape is recording

      /*******************************************************************************/
      /// @name Reversal

      void evaluate(); ///< Perform a full reverse evaluation of the tape.

      /*******************************************************************************/
      /// @name Reset

      void clearAdjoints();  ///< Clear all adjoint values and set them to zero.
      void reset(bool resetAdjoints = true);  ///< Reset the tape to the initial state for a fresh recording. See tutorial (TODO: ref) for remarks on repeated tape recording in CoDiPack.

      /*******************************************************************************/
      /// @name Tape information

      /// Default formatting of TapeValues
      template<typename Stream = std::ostream> void printStatistics(Stream& out = std::cout) const;
      /// Table header output of TapeValues
      template<typename Stream = std::ostream> void printTableHeader(Stream& out = std::cout) const;
      /// Table row output of TapeValues
      template<typename Stream = std::ostream> void printTableRow(Stream& out = std::cout) const;
      /// Get current tape values.
      TapeValues getTapeValues() const;
  };
}

