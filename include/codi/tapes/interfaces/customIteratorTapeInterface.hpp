/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2025 Chair for Scientific Computing (SciComp), University of Kaiserslautern-Landau
 * Homepage: http://scicomp.rptu.de
 * Contact:  Prof. Nicolas R. Gauger (codi@scicomp.uni-kl.de)
 *
 * Lead developers: Max Sagebaum, Johannes Blühdorn (SciComp, University of Kaiserslautern-Landau)
 *
 * This file is part of CoDiPack (http://scicomp.rptu.de/software/codi).
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
#include <memory>

#include "../../config.h"
#include "../../misc/macros.hpp"
#include "../data/position.hpp"
#include "../misc/lowLevelFunctionEntry.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief interface for callbacks from an custom tape iteration. See CustomIteratorTapeInterface for details.
   *
   * @tparam T_Real        The computation type of a tape, usually chosen as ActiveType::Real.
   * @tparam T_Identifier  The adjoint/tangent identification type of a tape, usually chosen as ActiveType::Identifier.
   */
  template<typename T_Real, typename T_Identifier>
  struct CallbacksInterface {
    public:
      using Real = CODI_DD(T_Real, double);           ///< See ReadWriteTapeInterface.
      using Identifier = CODI_DD(T_Identifier, int);  ///< See ReadWriteTapeInterface.

      using Tape = void;         ///< Any CoDiPack tape implementation.
      using EvalHandle = void*;  ///< See PrimalValueTapeTypes.

      /// Called for each statement in a Jacobian tape.
      void handleStatement(Identifier& lhsIndex, Config::ArgumentSize const& size, Real const* jacobians,
                           Identifier const* rhsIdentifiers);
      /// Called for each statement in a primal value tape.
      void handleStatement(EvalHandle const& evalHandle, Config::ArgumentSize const& nPassiveValues,
                           size_t& linearAdjointPosition, char* stmtData);
      /// Called for each low level function.
      void handleLowLevelFunction(LowLevelFunctionEntry<Tape, Real, Identifier> const& func, ByteDataView& llfData);
  };

  /**
   * @brief Iterate over the statement and low level function entries in a tape.
   *
   * Access to adjoint, primal and other tape data needs to be capture in the callbacks object. The callback object
   * needs to implement the CallbacksInterface.
   *
   * @tparam T_Position  Global tape position, usually chosen as Tape::Position.
   */
  template<typename T_Position>
  struct CustomIteratorTapeInterface {
    public:
      using Position = CODI_DD(T_Position, EmptyPosition);  ///< See CustomIteratorTapeInterface.

      /*******************************************************************************/
      /// @name Tape iteration and editing

      /// Iterate over the tape in a generalized fashion. callbacks needs to implement codi::CallbacksInterface.
      template<typename Callbacks>
      void iterateForward(Callbacks&& callbacks, Position const& start, Position const& end);

      /// Iterate over the tape in a generalized fashion. callbacks needs to implement codi::CallbacksInterface.
      template<typename Callbacks>
      void iterateForward(Callbacks&& callbacks);

      /// Iterate over the tape in a generalized fashion. callbacks needs to implement codi::CallbacksInterface.
      template<typename Callbacks>
      void iterateReverse(Callbacks&& callbacks, Position const& start, Position const& end);

      /// Iterate over the tape in a generalized fashion. callbacks needs to implement codi::CallbacksInterface.
      template<typename Callbacks>
      void iterateReverse(Callbacks&& callbacks);
  };
}
