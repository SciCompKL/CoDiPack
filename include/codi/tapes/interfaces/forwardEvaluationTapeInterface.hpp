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
#include "positionalEvaluationTapeInterface.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Forward AD evaluation of a recorded tape.
   *
   * See \ref TapeInterfaces for a general overview of the tape interface design in CoDiPack.
   *
   * Here is an example for a forward tape evaluation (documentation/examples/forwardEvaluationTapeInterface.cpp):
   * \snippet examples/forwardEvaluationTapeInterface.cpp Forward tape evaluation
   *
   * @tparam T_Position  Global tape position type, usually chosen as Tape::Position.
   */
  template<typename T_Position>
  struct ForwardEvaluationTapeInterface : public virtual PositionalEvaluationTapeInterface<T_Position> {
    public:

      using Position = CODI_DD(T_Position, EmptyPosition);  ///< See ForwardEvaluationTapeInterface.

      /*******************************************************************************/
      /// @name Interface definition

      /// Perform a forward evaluation of a part of the tape. It has to hold start <= end.
      void evaluateForward(Position const& start, Position const& end);

      /// Perform a forward evaluation of the full tape.
      void evaluateForward();
  };
}
