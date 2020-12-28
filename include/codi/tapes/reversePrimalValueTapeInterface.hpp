/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2021 Chair for Scientific Computing (SciComp), TU Kaiserslautern
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
 * Authors:
 *  - SciComp, TU Kaiserslautern:
 *     Max Sagebaum
 *     Tim Albring
 *     Johannes Blühdorn
 */


#pragma once

#include "reverseTapeInterface.hpp"

/**
 * @brief Global namespace for CoDiPack - Code Differentiation Package
 */
namespace codi {

  /**
   * @brief Interface common to all reverse tapes that use primal values.
   *
   * This extended interface each reverse tape with primal values has to implement. In addition to the
   * ReverseTapeInterface it provides function to manage the primal values of the tape.
   *
   * @tparam               Real  Floating point type of the gradients.
   * @tparam   GradientDataType  The data the tape uses to identify each active variable
   *                               and where the tape can store information about the
   *                               gradient.
   * @tparam   GradientValueType The value type that is used for the gradient calculation.
   * @tparam TapeImplementation  The implementing tape of the interface. It is needed to define the active type
   *                               for the registration of the variables.
   * @tparam           Position  Position used by the implementing tape.
   *
   */
  template <typename Real, typename GradientDataType, typename GradientValueType, typename TapeImplementation, typename Position>
  class ReversePrimalValueTapeInterface : virtual public ReverseTapeInterface<Real, GradientDataType, GradientValueType, TapeImplementation, Position> {
    public:

      /**
       * @brief Reverts the primal values to the given position.
       *
       * The method performs a reverse interpretation of the tape but only updates the primal
       * values with there overwritten values.
       *
       * Can be used to prepare a forward evaluation.
       *
       * @param[in] pos  The position to which the tape is reset.
       */
      virtual void revertPrimals(Position const& pos) = 0;

  };
}
