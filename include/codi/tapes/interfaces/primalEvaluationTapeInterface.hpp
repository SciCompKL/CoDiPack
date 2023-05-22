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
#include "positionalEvaluationTapeInterface.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Perform a primal reevaluation of the tape.
   *
   * See \ref TapeInterfaces for a general overview of the tape interface design in CoDiPack.
   *
   * Whether the tape manages primal values is indicated by the static constant `HasPrimalValues`.
   *
   * In a primal value tape, the correctness of the primal values is very important. The tapes should be programmed such
   * that the primal values stored in the tape are always up to date with the state of the program. Only through user
   * interaction this sync in states can be broken, but then the user should know what he is doing.
   *
   * The primal evaluation is used to reevaluate the primal values stored in the tape for different values of the
   * registered inputs. Note that this reevaluation follows the control flow that was observed during recording. The
   * control flow statements themselves, e.g. if constructs or loops, are not treated by CoDiPack. The user cannot
   * expect a reevaluation to choose different branches in if constructs or different numbers of loop iterations with
   * respect to the code that was recorded.
   *
   * Here is an example for a primal reevaluation (documentation/examples/primalEvaluationTapeInterface.cpp):
   * \snippet examples/primalEvaluationTapeInterface.cpp Primal evaluation
   *
   * @tparam T_Real        The computation type of a tape, usually chosen as ActiveType::Real.
   * @tparam T_Identifier  The adjoint/tangent identification of a tape, usually chosen as ActiveType::Identifier.
   * @tparam T_Position  Global tape position, usually chosen as Tape::Position.
   */
  template<typename T_Real, typename T_Identifier, typename T_Position>
  struct PrimalEvaluationTapeInterface : public virtual PositionalEvaluationTapeInterface<T_Position> {
    public:

      using Real = CODI_DD(T_Real, double);                 ///< See PrimalEvaluationTapeInterface.
      using Identifier = CODI_DD(T_Identifier, int);        ///< See PrimalEvaluationTapeInterface.
      using Position = CODI_DD(T_Position, EmptyPosition);  ///< See PrimalEvaluationTapeInterface.

      /*******************************************************************************/
      /// @name Interface definition

      static bool constexpr HasPrimalValues = CODI_UNDEFINED_VALUE;  ///< True if the tape has primal values.
      static bool constexpr RequiresPrimalRestore =
          CODI_UNDEFINED_VALUE;  ///< True if the primal state changes during a reverse or forward evaluation.

      /// Perform a partly (forward) reevaluation of the primals in the tape. It has to hold start <= end.
      void evaluatePrimal(Position const& start, Position const& end);

      /// Perform a full (forward) reevaluation of the primals in the tape.
      void evaluatePrimal();

      void setPrimal(Identifier const& identifier, Real const& gradient);  ///< Set primal value.
      Real const& getPrimal(Identifier const& identifier) const;           ///< Get primal value.

      Real& primal(Identifier const& identifier);              ///< Writable reference to primal value.
      Real const& primal(Identifier const& identifier) const;  ///< Read only reference to primal value.

      /// Revert the primals to the state indicated by pos.
      void revertPrimals(Position const& pos);
  };
}
