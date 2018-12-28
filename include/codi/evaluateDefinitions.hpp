/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2019 Chair for Scientific Computing (SciComp), TU Kaiserslautern
 * Homepage: http://www.scicomp.uni-kl.de
 * Contact:  Prof. Nicolas R. Gauger (codi@scicomp.uni-kl.de)
 *
 * Lead developers: Max Sagebaum, Tim Albring (SciComp, TU Kaiserslautern)
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
 * Authors: Max Sagebaum, Tim Albring, (SciComp, TU Kaiserslautern)
 */

#pragma once

#include "configure.h"
#include "typeTraits.hpp"
#include "tapeTypes.hpp"

/**
 * @brief Global namespace for CoDiPack - Code Differentiation Package
 */
namespace codi {

  /**
   * @brief Function definitions for primal value tapes.
   *
   * This are the functions that need to be implemented by primal value tapes, handle generators and expressions.
   *
   * @tparam ReverseTapeTypes  The basic type definitions for the tape. Need to define everything from ReverseTapeTypes.
   */
  template<typename ReverseTapeTypes>
  struct EvaluateDefinitions {

      CODI_INLINE_REVERSE_TAPE_TYPES(ReverseTapeTypes)

      /**
       * @brief Function definition for the primal evaluation of a statement in a tape.
       *
       * The function is called in the tapes for the primal evaluation of a statement. The arrays are still in the
       * global address space and need to be accessed with corresponding position.
       */
      typedef Real (*PrimalFunc)(const StatementInt& passiveActives,
                                 size_t& indexPos,
                                 Index* &indices,
                                 size_t& passivePos,
                                 Real* &passives,
                                 size_t& constantPos,
                                 PassiveReal* &constants,
                                 Real* primalVector);

      /**
       * @brief Function definition for the reverse evaluation of a statement in a tape.
       *
       * The function is called in the tapes for the reverse evaluation of a statement. The arrays are still in the
       * global address space and need to be accessed with corresponding position.
       */
      typedef void (*AdjointFunc)(const PRIMAL_SEED_TYPE& adj,
                                  const StatementInt& passiveActives,
                                  size_t& indexPos,
                                  Index* &indices,
                                  size_t& passivePos,
                                  Real* &passives,
                                  size_t& constantPos,
                                  PassiveReal* &constants,
                                  Real* primalVector,
                                  PRIMAL_ADJOINT_TYPE* adjoints);

      /**
       * @brief Function definition for the forward evaluation of a statement in a tape.
       *
       * The function is called in the tapes for the forward evaluation of a statement. The arrays are still in the
       * global address space and need to be accessed with corresponding position.
       */
      typedef Real (*TangentFunc)(const Real& adj,
                                  GradientValue& lhsAdj,
                                  const StatementInt& passiveActives,
                                  size_t& indexPos,
                                  Index* &indices,
                                  size_t& passivePos,
                                  Real* &passives,
                                  size_t& constantPos,
                                  PassiveReal* &constants,
                                  Real* primalVector,
                                  PRIMAL_ADJOINT_TYPE* adjoints);

      /**
       * @brief Function definition for the primal evaluation of a statement in an expression.
       *
       * The function is called in the expressions for the primal evaluation of a statement. The arrays are
       * defined in the local address space of the expression. They need to indexed starting with zero.
       */
      typedef Real (*PrimalExprFunc)(const Index* indices,
                                     const PassiveReal* constants,
                                     const Real* primalVector);

      /**
       * @brief Function definition for the reverse evaluation of a statement in an expression.
       *
       * The function is called in the expressions for the reverse evaluation of a statement. The arrays are
       * defined in the local address space of the expression. They need to indexed starting with zero.
       */
      typedef void (*AdjointExprFunc)(const Real& adj,
                                      const Index* indices,
                                      const PassiveReal* constants,
                                      const Real* primalVector,
                                      PRIMAL_ADJOINT_TYPE* adjoints);

      /**
       * @brief Function definition for the forward evaluation of a statement in an expression.
       *
       * The function is called in the expressions for the forward evaluation of a statement. The arrays are
       * defined in the local address space of the expression. They need to indexed starting with zero.
       */
      typedef Real (*TangentExprFunc)(const Real& adj,
                                      GradientValue& lhsAdj,
                                      const Index* indices,
                                      const PassiveReal* constants,
                                      const Real* primalVector,
                                      PRIMAL_ADJOINT_TYPE* adjoints);
  };
}
