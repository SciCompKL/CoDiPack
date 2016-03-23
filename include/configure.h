/*
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

#include "exceptions.hpp"
#include "macros.h"

/**
 * @brief Global namespace for CoDiPack - Code Differentiation Package
 */
namespace codi {

  /**
   * @brief Constant for the conversion from byte to megabyte.
   */
  const double BYTE_TO_MB = 1.0/1024.0/1024.0;

  /**
   * @brief Type for the maximum number of variables a operation can have.
   */
  typedef uint8_t StatementInt;

  #ifndef CODI_SmallChunkSize
    #define CODI_SmallChunkSize 32768
  #endif
  /**
   * @brief Default number of entries for all chunks that need smaller sizes.
   *
   * Default is 128 kb for 4 byte entries.
   */
  static size_t DefaultSmallChunkSize = CODI_SmallChunkSize; //TODO: Find optimal value
  #undef CODI_SmallChunkSize

  #ifndef CODI_ChunkSize
    #define CODI_ChunkSize 2097152
  #endif
  /**
   * @brief Default number of entries for all chunks.
   *
   * Default is 24 Mb for 12 byte entries.
   */
  static size_t DefaultChunkSize = CODI_ChunkSize; //TODO: Find optimal value
  #undef CODI_ChunkSize

  #ifndef CODI_UseMemsetInChunks
    #define CODI_UseMemsetInChunks true
  #endif
  /**
   * @brief Enables the use of memset to preallocate the memory.
   *
   * Modern systems can initialize the memory acquired with malloc or calloc
   * in a lazy fashion. The memory is then acquired on first use which can cause
   * performance issues. If the data is set to zero with memset the system will
   * directly allocate all the memory.
   */
  const bool UseMemsetInChunks = CODI_UseMemsetInChunks;
  #undef CODI_UseMemsetInChunks

  #ifndef CODI_CheckExpressionArguments
    #define CODI_CheckExpressionArguments false
  #endif
  /**
   * @brief Check if the arguments are inside the differentiable domain.
   *
   * The check enables for all function the validation of the arguments for
   * gradient evaluation. If the arguments are not valid a CoDiPack exceptions is
   * generated.
   */
  const bool CheckExpressionArguments = CODI_CheckExpressionArguments;
  #undef CODI_CheckExpressionArguments


  #ifndef CODI_OptIgnoreInvalidJacobies
    #define CODI_OptIgnoreInvalidJacobies false
  #endif
  /**
   * @brief Tapes push jacobies only if they are valid values.
   *
   * The check is used in the tape 'pushJacobi' function to disable the pushing of the
   * jacobies if they are nan or inf.
   */
  const bool OptIgnoreInvalidJacobies = CODI_OptIgnoreInvalidJacobies;
  #undef CODI_OptIgnoreInvalidJacobies


  #ifndef CODI_OptJacobiIsZero
    #define CODI_OptJacobiIsZero true
  #endif
  /**
   * @brief Tapes push jacobies only if they are none zero.
   *
   * The check is used in the tape 'pushJacobi' function to disable the pushing of the
   * jacobies if they are zero.
   */
  const bool OptJacobiIsZero = CODI_OptJacobiIsZero;
  #undef CODI_OptJacobiIsZero

  #ifndef CODI_OptCheckZeroIndex
    #define CODI_OptCheckZeroIndex true
  #endif
  /**
   * @brief Tapes push jacobies only if there index is not zero
   *
   * The check is used in the tape 'pushJacobi' and 'store' functions to disable the pushing of the
   * jacobies if there index is zero.
   */
  const bool OptCheckZeroIndex = CODI_OptCheckZeroIndex;
  #undef CODI_OptCheckZeroIndex

  #ifndef CODI_OptCheckEmptyStatements
    #define CODI_OptCheckEmptyStatements true
  #endif
  /**
   * @brief Tapes push statements only if at least one jacobi was pushed.
   *
   * The check is used in the tape 'store' function to disable the pushing of the
   * statement if no jacobi was pushed.
   */
  const bool OptCheckEmptyStatements = CODI_OptCheckEmptyStatements;
  #undef CODI_OptCheckEmptyStatements

  #ifndef CODI_OptTapeActivity
    #define CODI_OptTapeActivity true
  #endif
  /**
   * @brief Tapes can be disable for regions which do not need to be taped.
   *
   * If the option is set to true a tape can be enable or disabled which can be used to
   * disable the tape for code parts which do not need to be taped. If the option is set
   * to false the tape will always be active.
   */
  const bool OptTapeActivity = CODI_OptTapeActivity;
  #undef CODI_OptTapeActivity

  #ifndef CODI_OptZeroAdjoint
    #define CODI_OptZeroAdjoint true
  #endif
  /**
   * @brief Omits the evaluation of jacobies which are zero in the reverse sweep.
   *
   * If an adjoint seed is zero during the reverse sweep, all the updates for the
   * adjoint vector will be zero. Therefore the loop does not need to be evaluated.
   */
  const bool OptZeroAdjoint = CODI_OptZeroAdjoint;
  #undef CODI_OptZeroAdjoint

  #ifndef CODI_DisableAssignOptimization
    #define CODI_DisableAssignOptimization false
  #endif
  /**
   * @brief Disables the assign optimization for linear index tapes.
   *
   * An assign statement usually does not need to written for tapes
   * that use a linear increasing index scheme. The correspoinding
   * entry on the tape would just add the accumulated values for
   * the lhs to the rhs. This optimization can be dissabled with
   * this switch.
   */
  const bool OptDisableAssignOptimization = CODI_DisableAssignOptimization;
  #undef CODI_DisableAssignOptimization

  #ifndef CODI_AdjointHandle
    #define CODI_AdjointHandle false
  #endif
  #if CODI_AdjointHandle
    /**
     * @brief A function that is called for every statement that is written on the
     *        Jacobie tapes.
     *
     * The function can be used to extract information from the taping process.
     *
     * @param[in]      value  The primal value of the statement.
     * @param[in]   lhsIndex  The index on the left hand side of the statement, that
     *                        will be recorded.
     * @param[in]   jacobies  The pointer to the array of the Jacobies that will be stored
     *                        for the statement. jacobies[0] is the first argument, jacobies[1]
     *                        the second, etc.
     * @param[in] rhsIndices  The pointer to the array of the indices that will be stored
     *                        for the statement. rhsIndices[0] is the first argument, rhsIndices[1]
     *                        the second, etc.
     * @param[in]       size  The number of arguments that are stored for the statement.
     *
     * @tparam      Real  The type of the floating point values that is used in the tape.
     * @tparam IndexType  The type of the indices that is used in the tape.
     */
    template<typename Real, typename IndexType>
    void handleAdjointOperation(const Real& value, const IndexType lhsIndex, const Real* jacobies, const IndexType* rhsIndices, const int size);
  #endif

  #ifndef CODI_EnableAssert
    #define CODI_EnableAssert false
  #endif
  #ifndef codiAssert
    #if CODI_EnableAssert
      /**
       * @brief The assert function for CoDiPack it can be enabled with the preprocessor macro CODI_EnableAssert=true
       *
       * @param x The expression that is checked in the assert.
       */
      #define codiAssert(x) codi::checkAndOutputAssert(x, CODI_TO_STRING(x), __PRETTY_FUNCTION__, __FILE__, __LINE__)
    #else
      /**
       * @brief The assert function for CoDiPack it can be enabled with the preprocessor macro CODI_EnableAssert=true
       *
       * @param x The expression that is checked in the assert.
       */
      #define codiAssert(x) /* disabled by CODI_EnableAssert */
    #endif
  #endif

  /**
   * @brief Enable the check only if the option is set
   *
   * The macro ca be used to surround a code block with an if statement. If the option is set to true
   * the condition is evaluated and only if the condition is true the block after the macro is
   * executed. If the option is false the block will always be executed.
   *
   * @param    option  A constant global boolean. Only than the compiler can optimize the statement.
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
