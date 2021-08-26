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
#include "internalStatementRecordingTapeInterface.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Minimum tape interface for a working reverse tape implementation.
   *
   * See \ref TapeInterfaces for a general overview of the tape interface design in CoDiPack.
   *
   * For an example on how to uses this interface to perform an AD reverse mode recording and evaluation of a program,
   * please see tutorials \ref Tutorial_02_Reverse_mode_AD and \ref Tutorial_05_Repeated_tape_recordings).
   *
   * Implementation hints:
   * A tape should only record information if it is active, that is, everything between a call to
   * setActive() and setPassive(). A call to setActive() does not reset the tape in CoDiPack. A reset can only be
   * performed by a call to reset(). Hence, the user may skip unnecessarys parts of the recording by setting the tape
   * passive for these regions.
   *
   * Here is an example for using a tape (documentation/examples/reverseModeAD.cpp):
   * \snippet examples/reverseModeAD.cpp Reverse mode AD
   *
   * @tparam T_Real        The computation type of a tape, usually chosen as ActiveType::Real.
   * @tparam T_Gradient    The gradient type of a tape, usually chosen as ActiveType::Gradient.
   * @tparam T_Identifier  The adjoint/tangent identification of a tape, usually chosen as ActiveType::Identifier.
   */
  template<typename T_Real, typename T_Gradient, typename T_Identifier>
  struct ReverseTapeInterface : public virtual InternalStatementRecordingTapeInterface<T_Identifier>,
                                public virtual GradientAccessTapeInterface<T_Gradient, T_Gradient> {
    public:
      using Real = CODI_DD(T_Real, double);           ///< See ReverseTapeInterface.
      using Gradient = CODI_DD(T_Gradient, double);   ///< See ReverseTapeInterface.
      using Identifier = CODI_DD(T_Identifier, int);  ///< See ReverseTapeInterface.

      using PassiveReal = RealTraits::PassiveReal<Real>;  ///< Basic computation type.

      /*******************************************************************************/
      /// @name Recording

      /// Mark a value as input (independent) and make it active.
      template<typename Lhs>
      void registerInput(LhsExpressionInterface<Real, Gradient, ReverseTapeInterface, Lhs>& value);
      /// Mark a value as output (dependent).
      template<typename Lhs>
      void registerOutput(LhsExpressionInterface<Real, Gradient, ReverseTapeInterface, Lhs>& value);

      void setActive();       ///< Start/continue recording of statements.
      void setPassive();      ///< Stop/interrupt recording of statements.
      bool isActive() const;  ///< Check if the tape is recording.

      /*******************************************************************************/
      /// @name Reversal

      void evaluate();  ///< Perform a full reverse evaluation of the tape.

      /*******************************************************************************/
      /// @name Reset

      void clearAdjoints();                   ///< Clear all adjoint values, that is, set them to zero.
      void reset(bool resetAdjoints = true);  ///< Reset the tape to the initial state for a fresh recording. See
                                              ///< \ref Tutorial_05_Repeated_tape_recordings for remarks on repeated
                                              ///< tape recording in CoDiPack.

      /*******************************************************************************/
      /// @name Tape information

      /// Default formatting of TapeValues.
      template<typename Stream = std::ostream>
      void printStatistics(Stream& out = std::cout) const;
      /// Table header output of TapeValues.
      template<typename Stream = std::ostream>
      void printTableHeader(Stream& out = std::cout) const;
      /// Table row output of TapeValues.
      template<typename Stream = std::ostream>
      void printTableRow(Stream& out = std::cout) const;
      /// Get current tape values.
      TapeValues getTapeValues() const;
  };
}
