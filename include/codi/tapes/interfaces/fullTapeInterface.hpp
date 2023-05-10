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
#include "customAdjointVectorEvaluationTapeInterface.hpp"
#include "dataManagementTapeInterface.hpp"
#include "externalFunctionTapeInterface.hpp"
#include "forwardEvaluationTapeInterface.hpp"
#include "gradientAccessTapeInterface.hpp"
#include "identifierInformationTapeInterface.hpp"
#include "internalStatementRecordingTapeInterface.hpp"
#include "manualStatementPushTapeInterface.hpp"
#include "positionalEvaluationTapeInterface.hpp"
#include "preaccumulationEvaluationTapeInterface.hpp"
#include "primalEvaluationTapeInterface.hpp"
#include "reverseTapeInterface.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Full tape interface that supports all features of CoDiPack.
   *
   * See \ref TapeInterfaces for a general overview of the tape interface design in CoDiPack.
   *
   * A tape that implements this interface correctly can be used in all helper structures of CoDiPack.
   *
   * @tparam T_Real        The computation type of a tape, usually chosen as ActiveType::Real.
   * @tparam T_Gradient    The gradient type of a tape, usually chosen as ActiveType::Gradient.
   * @tparam T_Identifier  The adjoint/tangent identification type of a tape, usually chosen as ActiveType::Identifier.
   * @tparam T_Position  Global tape position, usually chosen as Tape::Position.
   */
  template<typename T_Real, typename T_Gradient, typename T_Identifier, typename T_Position>
  struct FullTapeInterface
      : public virtual CustomAdjointVectorEvaluationTapeInterface<T_Position>,
        public virtual DataManagementTapeInterface<T_Real, T_Identifier>,
        public virtual ExternalFunctionTapeInterface<T_Real, T_Gradient, T_Identifier>,
        public virtual ForwardEvaluationTapeInterface<T_Position>,
        public virtual GradientAccessTapeInterface<T_Gradient, T_Identifier>,
        public virtual IdentifierInformationTapeInterface<T_Real, T_Gradient, T_Identifier>,
        public virtual InternalStatementRecordingTapeInterface<T_Identifier>,
        public virtual ManualStatementPushTapeInterface<T_Real, T_Gradient, T_Identifier>,
        public virtual PositionalEvaluationTapeInterface<T_Position>,
        public virtual PreaccumulationEvaluationTapeInterface<T_Real, T_Gradient, T_Identifier, T_Position>,
        public virtual PrimalEvaluationTapeInterface<T_Real, T_Identifier, T_Position>,
        public virtual ReverseTapeInterface<T_Real, T_Gradient, T_Identifier> {
    public:

      using Real = CODI_DD(T_Real, double);                 ///< See FullTapeInterface.
      using Gradient = CODI_DD(T_Gradient, double);         ///< See FullTapeInterface.
      using Identifier = CODI_DD(T_Identifier, int);        ///< See FullTapeInterface.
      using Position = CODI_DD(T_Position, EmptyPosition);  ///< See FullTapeInterface.

#if CODI_IDE
      using TapeTypes = CODI_ANY;
#endif

      using PositionalEvaluationTapeInterface<T_Position>::clearAdjoints;
      using ReverseTapeInterface<T_Real, T_Gradient, T_Identifier>::clearAdjoints;

      using PositionalEvaluationTapeInterface<T_Position>::evaluate;
      using ReverseTapeInterface<T_Real, T_Gradient, T_Identifier>::evaluate;
      using CustomAdjointVectorEvaluationTapeInterface<T_Position>::evaluate;

      using ForwardEvaluationTapeInterface<T_Position>::evaluateForward;
      using CustomAdjointVectorEvaluationTapeInterface<T_Position>::evaluateForward;
  };
}
