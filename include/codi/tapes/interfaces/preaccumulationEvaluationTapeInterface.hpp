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
#include "forwardEvaluationTapeInterface.hpp"
#include "manualStatementPushTapeInterface.hpp"
#include "positionalEvaluationTapeInterface.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Perform tape evaluations but ensure that the state prior to evaluation equals the state after evaluation.
   *
   * See \ref TapeInterfaces for a general overview of the tape interface design in CoDiPack.
   *
   * These interface functions are used for small tape evaluation where only a part of the tape is evaluated. Especially
   * for primal value tapes, it is essential that the primal value store in the background is in sync with the program
   * state. The normal evaluate methods in these tapes copy the primal value vector and perform all operations on the
   * copied vector. In the *KeepState() methods, they ensure the correctness of the primal value vector by performing
   * e.g. a primal evaluation after the reverse evaluation. This will yield better performance for small tape ranges in
   * the evaluation.
   *
   * @tparam T_Real        The computation type of a tape, usually chosen as ActiveType::Real.
   * @tparam T_Gradient    The gradient type of a tape usually, chosen as ActiveType::Gradient.
   * @tparam T_Identifier  The adjoint/tangent identification type of a tape, usually chosen as ActiveType::Identifier.
   * @tparam T_Position  Global tape position, usually chosen as Tape::Position.
   */
  template<typename T_Real, typename T_Gradient, typename T_Identifier, typename T_Position>
  struct PreaccumulationEvaluationTapeInterface
      : public virtual PositionalEvaluationTapeInterface<T_Position>,
        public virtual ForwardEvaluationTapeInterface<T_Position>,
        public virtual ManualStatementPushTapeInterface<T_Real, T_Gradient, T_Identifier> {
    public:

      using Real = CODI_DD(T_Real, double);                 ///< See PreaccumulationEvaluationTapeInterface.
      using Gradient = CODI_DD(T_Gradient, double);         ///< See PreaccumulationEvaluationTapeInterface.
      using Identifier = CODI_DD(T_Identifier, int);        ///< See PreaccumulationEvaluationTapeInterface.
      using Position = CODI_DD(T_Position, EmptyPosition);  ///< See PreaccumulationEvaluationTapeInterface.

      /*******************************************************************************/
      /// @name Interface definition

      /**
       * @brief Perform a tape evaluation but restore the state afterwards such that it is the same as when the
       * evaluation started. It hast to hold start >= end.
       *
       * Automatic adjoints management involves bounds checking, resizing, and locking, see AdjointsManagement for
       * details.
       */
      void evaluateKeepState(Position const& start, Position const& end,
                             AdjointsManagement adjointsManagement = AdjointsManagement::Automatic);

      /**
       * @brief Perform a tape evaluation but restore the state afterwards such that it is the same as when the
       * evaluation started. It hast to hold start <= end.
       *
       * Automatic adjoints management involves bounds checking, resizing, and locking, see AdjointsManagement for
       * details.
       */
      void evaluateForwardKeepState(Position const& start, Position const& end,
                                    AdjointsManagement adjointsManagement = AdjointsManagement::Automatic);
  };
}
