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

#include "../../config.h"
#include "../../misc/macros.hpp"
#include "../data/position.hpp"
#include "../misc/tapeParameters.hpp"

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

      /**
       * @brief Perform a reverse evaluation for a part of the tape. It hast to hold start >= end.
       *
       * Automatic adjoints management involves bounds checking, resizing, and locking, see AdjointsManagement for
       * details.
       */
      void evaluate(Position const& start, Position const& end,
                    AdjointsManagement adjointsManagement = AdjointsManagement::Automatic);

      /**
       * @brief Clear all adjoints that would be set in a tape evaluation from start to end. It has to hold start >=
       * end.
       *
       * Automatic adjoints management involves locking, see AdjointsManagement for details.
       */
      void clearAdjoints(Position const& start, Position const& end,
                         AdjointsManagement adjointsManagement = AdjointsManagement::Automatic);

      Position getPosition() const;      ///< Current position of the tape.
      Position getZeroPosition() const;  ///< Initial position of the tape.

      /**
       * @brief Reset the tape to the provided position.
       *
       * Automatic adjoints management involves locking, see AdjointsManagement for details.
       */
      void resetTo(Position const& pos, bool resetAdjoints = true,
                   AdjointsManagement adjointsManagement = AdjointsManagement::Automatic);
  };
}
