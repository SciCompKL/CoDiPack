/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2023 Chair for Scientific Computing (SciComp), University of Kaiserslautern-Landau
 * Homepage: http://www.scicomp.uni-kl.de
 * Contact:  Prof. Nicolas R. Gauger (codi@scicomp.uni-kl.de)
 *
 * Lead developers: Max Sagebaum, Johannes Blühdorn (SciComp, University of Kaiserslautern-Landau)
 *
 * This file is part of CoDiPack (http://www.scicomp.uni-kl.de/software/codi).
 *
 * CoDiPack is free software: you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * CoDiPack is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty
 * of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 * See the GNU General Public License for more details.
 * You should have received a copy of the GNU
 * General Public License along with CoDiPack.
 * If not, see <http://www.gnu.org/licenses/>.
 *
 * For other licensing options please contact us.
 *
 * Authors:
 *  - SciComp, University of Kaiserslautern-Landau:
 *    - Max Sagebaum
 *    - Johannes Blühdorn
 *    - Former members:
 *      - Tim Albring
 */
#pragma once

#include <iostream>

#include "../../config.h"
#include "../../expressions/logic/traversalLogic.hpp"
#include "../../misc/macros.hpp"
#include "../../traits/expressionTraits.hpp"
#include "../../traits/realTraits.hpp"
#include "../misc/tapeParameters.hpp"
#include "../misc/tapeValues.hpp"
#include "gradientAccessTapeInterface.hpp"
#include "internalStatementRecordingTapeInterface.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  template<typename T_Real, typename T_Gradient, typename T_Tape, typename T_Impl>
  struct LhsExpressionInterface;

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
                                public virtual GradientAccessTapeInterface<T_Gradient, T_Identifier> {
    public:
      using Real = CODI_DD(T_Real, double);           ///< See ReverseTapeInterface.
      using Gradient = CODI_DD(T_Gradient, double);   ///< See ReverseTapeInterface.
      using Identifier = CODI_DD(T_Identifier, int);  ///< See ReverseTapeInterface.

      using PassiveReal = RealTraits::PassiveReal<Real>;  ///< Basic computation type.

      /*******************************************************************************/
      /// @name Recording

      /// Mark a value as input (independent) and make it active.
      ///
      /// @tparam Lhs  Class that implements the LhsExpressionInterface. See also LhsExpressionInterface.
      /// @tparam Tape  Tape implementation used in the LhsExpressionInterface. See also LhsExpressionInterface.
      template<typename Lhs, typename Tape>
      void registerInput(LhsExpressionInterface<Real, Gradient, Tape, Lhs>& value);
      /// Mark a value as output (dependent).
      ///
      /// @tparam Lhs  Class that implements the LhsExpressionInterface. See also LhsExpressionInterface.
      /// @tparam Tape  Tape implementation used in the LhsExpressionInterface. See also LhsExpressionInterface.
      template<typename Lhs, typename Tape>
      void registerOutput(LhsExpressionInterface<Real, Gradient, Tape, Lhs>& value);

      void setActive();       ///< Start/continue recording of statements.
      void setPassive();      ///< Stop/interrupt recording of statements.
      bool isActive() const;  ///< Check if the tape is recording.

      /*******************************************************************************/
      /// @name Reversal

      /**
       * @brief Perform a full reverse evaluation of the tape.
       *
       * Automatic adjoints management involves bounds checking, resizing, and locking, see AdjointsManagement for
       * details.
       */
      void evaluate(AdjointsManagement adjointsManagement = AdjointsManagement::Automatic);

      /*******************************************************************************/
      /// @name Reset

      /**
       * @brief Clear all adjoint values, that is, set them to zero.
       *
       * Automatic adjoints management involves locking, see AdjointsManagement for details.
       */
      void clearAdjoints(AdjointsManagement adjointsManagement = AdjointsManagement::Automatic);

      /**
       * @brief Reset the tape to the initial state for a fresh recording.
       *
       * See \ref Tutorial_05_Repeated_tape_recordings for remarks on repeated tape recording in CoDiPack.
       *
       * Automatic adjoints management involves locking, see AdjointsManagement for details.
       */
      void reset(bool resetAdjoints = true, AdjointsManagement adjointsManagement = AdjointsManagement::Automatic);

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
