/**
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015 Chair for Scientific Computing, TU Kaiserslautern
 *
 * This file is part of CoDiPack.
 *
 * CoDiPack is free software: you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation, either version 2 of the
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
 * Authors: TODO
 */
#pragma once

#include <stdint.h>


#define NDEBUG
#include <assert.h>

namespace codi {

  /**
   * @brief Type for the maximum number of variables a opeartion can have.
   */
  typedef uint8_t OperationInt;

  /**
   * @brief Default size for all chunks.
   */
  static size_t DefaultChunkSize = 13107200; //TODO: Find optimal value


  /**
   * @brief Tapes push jacobies only if they are none zero.
   *
   * The check is used in the tape 'pushJacobi' function to disable the pushing of the
   * jacobies if they are zero.
   */
  const bool OptJacobiIsZero = true;

  /**
   * @brief Tapes can be disable for regions which do not need to be taped.
   *
   * If the option is set to true a tape can be enable or disabled which can be used to
   * disable the tape for code parts which do not need to be taped. If the option is set
   * to false the tape will always be active.
   */
  const bool OptTapeActivity = true;

  /**
   * @brief Omits the evaluation of jacobies which are zero in the reverse sweep.
   *
   * If an adjoint seed is zero during the reverse sweep, all the updates for the
   * adjoint vector will be zero. Thefore the loop does not need to be evaluted.
   */
  const bool OptZeroAdjoint = true;

  /**
   * @brief Enable the check only if the option is set
   *
   * The macro ca be used to surrund a code block with an if statement. If the option is set to true
   * the condition is evaluated and only if the condition is true the block after the macro is
   * executed. If the option is false the block will always be executed.
   *
   * @param option     A constant global boolean. Only than the compiler can optimize the statement.
   * @param condition  The condition which is only evaluted if 'option' is set to true.
   */
#  define ENABLE_CHECK(option, condition) if(!(option) || (condition))

}
