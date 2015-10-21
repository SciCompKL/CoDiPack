/**
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015 Chair for Scientific Computing (SciComp), TU Kaiserslautern
 * Homepage: http://www.scicomp.uni-kl.de
 * Contact:  Prof. Nicolas R. Gauger (codi@scicomp.uni-kl.de)
 *
 * Lead developers: Max Sagebaum, Tim Albring (SciComp, TU Kaiserslautern)
 *
 * This file is part of CoDiPack (http://www.scicomp.uni-kl.de/software/codi).
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
 * Authors: Max Sagebaum, Tim Albring, (SciComp, TU Kaiserslautern)
 */

#pragma once

#include <stdint.h>


#define NDEBUG
#include <assert.h>

namespace codi {

  /**
   * @brief Type for the maximum number of variables a operation can have.
   */
  typedef uint8_t StatementInt;

  /**
   * @brief Default number of entries for all chunks.
   */
  #ifndef CODI_ChunkSize
    #define CODI_ChunkSize 13107200
  #endif
  static size_t DefaultChunkSize = CODI_ChunkSize; //TODO: Find optimal value

  /**
   * @brief Enables the use of memset to preallocate the memory.
   *
   * Modern systems can initalize the memory acquired with malloc or calloc
   * in a lazy fashion. The memory is then acquired on first use which can cause
   * performace issues. If the data is set to zero with memset the system will
   * directly allocate all the memory.
   */
  #ifndef CODI_UseMemsetInChunks
    #define CODI_UseMemsetInChunks true
  #endif
  const bool UseMemsetInChunks = CODI_UseMemsetInChunks;
  #undef CODI_UseMemsetInChunks

  /**
   * @brief Check if the arguments are inside the differentiable domain.
   *
   * The check enables for all function the validation of the arguments for
   * gradient evaluation. If the arguments are not valid a CoDiPack exceptions is
   * generated.
   */
  #ifndef CODI_CheckExpressionArguments
    #define CODI_CheckExpressionArguments false
  #endif
  const bool CheckExpressionArguments = CODI_CheckExpressionArguments;
  #undef CODI_CheckExpressionArguments


  /**
   * @brief Tapes push jacobies only if they are valid values.
   *
   * The check is used in the tape 'pushJacobi' function to disable the pushing of the
   * jacobies if they are nan or inf.
   */
  #ifndef CODI_OptIgnoreInvalidJacobies
    #define CODI_OptIgnoreInvalidJacobies false
  #endif
  const bool OptIgnoreInvalidJacobies = CODI_OptIgnoreInvalidJacobies;
  #undef CODI_OptIgnoreInvalidJacobies


  /**
   * @brief Tapes push jacobies only if they are none zero.
   *
   * The check is used in the tape 'pushJacobi' function to disable the pushing of the
   * jacobies if they are zero.
   */
  #ifndef CODI_OptJacobiIsZero
    #define CODI_OptJacobiIsZero true
  #endif
  const bool OptJacobiIsZero = CODI_OptJacobiIsZero;
  #undef CODI_OptJacobiIsZero

  /**
   * @brief Tapes can be disable for regions which do not need to be taped.
   *
   * If the option is set to true a tape can be enable or disabled which can be used to
   * disable the tape for code parts which do not need to be taped. If the option is set
   * to false the tape will always be active.
   */
  #ifndef CODI_OptTapeActivity
    #define CODI_OptTapeActivity true
  #endif
  const bool OptTapeActivity = CODI_OptTapeActivity;
  #undef CODI_OptTapeActivity

  /**
   * @brief Omits the evaluation of jacobies which are zero in the reverse sweep.
   *
   * If an adjoint seed is zero during the reverse sweep, all the updates for the
   * adjoint vector will be zero. Therefore the loop does not need to be evaluated.
   */
  #ifndef CODI_OptZeroAdjoint
    #define CODI_OptZeroAdjoint true
  #endif
  const bool OptZeroAdjoint = CODI_OptZeroAdjoint;
  #undef CODI_OptZeroAdjoint

  /**
   * @brief Enable the check only if the option is set
   *
   * The macro ca be used to surround a code block with an if statement. If the option is set to true
   * the condition is evaluated and only if the condition is true the block after the macro is
   * executed. If the option is false the block will always be executed.
   *
   * @param option     A constant global boolean. Only than the compiler can optimize the statement.
   * @param condition  The condition which is only evaluated if 'option' is set to true.
   */
#  define ENABLE_CHECK(option, condition) if(!(option) || (condition))

  /**
   * @brief Needed to disable warnings about unused parameters.
   *
   * Is also necessary because of doxygen parameter handling.
   */
#  define CODI_UNUSED(name) (void)name
}
