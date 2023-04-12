/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2022 Chair for Scientific Computing (SciComp), TU Kaiserslautern
 * Homepage: http://www.scicomp.uni-kl.de
 * Contact:  Prof. Nicolas R. Gauger (codi@scicomp.uni-kl.de)
 *
 * Lead developers: Max Sagebaum, Johannes Blühdorn (SciComp, TU Kaiserslautern)
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
 *  - SciComp, TU Kaiserslautern:
 *    - Max Sagebaum
 *    - Johannes Blühdorn
 *    - Former members:
 *      - Tim Albring
 */
#pragma once

#include "../../misc/macros.hpp"
#include "../../config.h"
#include "../data/position.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Reverse AD evaluation for parts of a recorded tape.
   *
   * See \ref TapeInterfaces for a general overview of the tape interface design in CoDiPack.
   *
   * The definitions in the ReverseTapeInterface provide only methods that operate on the full tape. With the methods
   * in this interface, all these operations can be performed on parts of the tape.
   *
   * Here is an example for a positional tape evaluation (\ref Example_03_Positional_tape_evaluations):
   * \snippet examples/Example_03_Positional_tape_evaluations.cpp Positional evaluation
   *
   * @tparam T_Position  Global tape position, usually chosen as Tape::Position.
   */
  template<typename T_Position>
  struct PositionalEvaluationTapeInterface {
    public:

      using Position = CODI_DD(T_Position, EmptyPosition);  ///< See PositionalEvaluationTapeInterface.

      /*******************************************************************************/
      /// @name Interface definition

      /// Perform a reverse evaluation for a part of the tape. It hast to hold start >= end.
      void evaluate(Position const& start, Position const& end);

      /// Clear all adjoints that would be set in a tape evaluation from start to end. It has to hold start >= end.
      void clearAdjoints(Position const& start, Position const& end);

      Position getPosition() const;      ///< Current position of the tape.
      Position getZeroPosition() const;  ///< Initial position of the tape.

      void resetTo(Position const& pos, bool resetAdjoints = true);  ///< Rest the tape to the provided position.
  };
}
