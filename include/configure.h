/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2017 Chair for Scientific Computing (SciComp), TU Kaiserslautern
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

#include <stdint.h>

#include "exceptions.hpp"
#include "macros.h"

/**
 * @brief Global namespace for CoDiPack - Code Differentiation Package
 */
namespace codi {

  #define CODI_MAJOR_VERSION 1
  #define CODI_MINOR_VERSION 4
  #define CODI_BUILD_VERSION 1
  #define CODI_VERSION "1.4.1"

  /**
   * @brief Constant for the conversion from byte to megabyte.
   */
  const double BYTE_TO_MB = 1.0/1024.0/1024.0;

   /**
   * @brief Macro for forcing the inlining of the expression templates.
   *
   * The macro defines the attribute of the function such that it is directly inlined
   * and not just an recommendation for the compiler.
   *
   * Currently it is defined for intel and gcc.
   */
  #ifndef CODI_UseForcedInlines
    #define CODI_UseForcedInlines 0
  #endif
  #if CODI_UseForcedInlines
    #if defined(__INTEL_COMPILER) | defined(_MSC_VER)
      #define CODI_INLINE __forceinline
    #elif defined(__GNUC__)
      #define CODI_INLINE inline __attribute__((always_inline))
    #else
      #warning Could not determine compiler for forced inline definitions. Using inline.
      #define CODI_INLINE inline
    #endif
  #else
    #define CODI_INLINE inline
  #endif
  #undef CODI_UseForcedInlines

   /**
   * @brief Macro for avoiding the inlining of function.
   *
   * The macro defines the attribute of the function such that it is no longer considered for inlining.
   *
   * It is defined as an function attribute.
   */
  #ifndef CODI_UseAvoidedInlines
    #define CODI_UseAvoidedInlines 1
  #endif
  #if CODI_UseAvoidedInlines
    #if defined(_MSC_VER)
      #define CODI_NO_INLINE __declspec(noinline)
    #else
      #define CODI_NO_INLINE __attribute__((noinline))
    #endif
  #else
    #define CODI_NO_INLINE /* no avoiding of inline defined */
  #endif
  #undef CODI_UseAvoidedInlines


  /**
   * @brief Type for the maximum number of variables a operation can have.
   */
  typedef uint8_t StatementInt;

  /**
   * @brief The maximum value of a statement int.
   */
  const size_t MaxStatementIntSize = 256;


  #ifndef CODI_SmallChunkSize
    #define CODI_SmallChunkSize 32768
  #endif
  /**
   * @brief Default number of entries for all chunks that need smaller sizes.
   *
   * Default is 128 kb for 4 byte entries.
   *
   * It can be set with the preprocessor macro CODI_SmallChunkSize=<size>
   */
  static size_t DefaultSmallChunkSize = CODI_SmallChunkSize;
  #undef CODI_SmallChunkSize

  #ifndef CODI_ChunkSize
    #define CODI_ChunkSize 2097152
  #endif
  /**
   * @brief Default number of entries for all chunks.
   *
   * Default is 24 Mb for 12 byte entries.
   *
   * It can be set with the preprocessor macro CODI_ChunkSize=<size>
   */
  static size_t DefaultChunkSize = CODI_ChunkSize;
  #undef CODI_ChunkSize

  #ifndef CODI_CheckExpressionArguments
    #define CODI_CheckExpressionArguments false
  #endif
  /**
   * @brief Check if the arguments are inside the differentiable domain.
   *
   * The check enables for all function the validation of the arguments for
   * gradient evaluation. If the arguments are not valid a CoDiPack exceptions is
   * generated.
   *
   * It can be set with the preprocessor macro CODI_CheckExpressionArguments=<true/false>
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
   *
   * It can be set with the preprocessor macro CODI_OptIgnoreInvalidJacobies=<true/false>
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
   *
   * It can be set with the preprocessor macro CODI_OptJacobiIsZero=<true/false>
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
   *
   * It can be set with the preprocessor macro CODI_OptCheckZeroIndex=<true/false>
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
   *
   * It can be set with the preprocessor macro CODI_OptCheckEmptyStatements=<true/false>
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
   *
   * It can be set with the preprocessor macro CODI_OptTapeActivity=<true/false>
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
   *
   * It can be set with the preprocessor macro CODI_OptZeroAdjoint=<true/false>
   */
  const bool OptZeroAdjoint = CODI_OptZeroAdjoint;
  #undef CODI_OptZeroAdjoint

  #ifndef CODI_DisableAssignOptimization
    #define CODI_DisableAssignOptimization false
  #endif
  /**
   * @brief Disables the assign optimization for linear index tapes.
   *
   * An assign statement usually does not need to be written for tapes
   * that use a linear increasing index scheme. The corresponding
   * entry on the tape would just add the accumulated values for
   * the lhs to the rhs. This optimization can be disabled with
   * this switch.
   *
   * It can be set with the preprocessor macro CODI_DisableAssignOptimization=<true/false>
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
     * It can be set with the preprocessor macro CODI_AdjointHandle=<true/false>
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
     * @tparam      Real  The type of the floating point values that are used in the tape.
     * @tparam IndexType  The type of the indices that are used in the tape.
     */
    template<typename Real, typename IndexType>
    void handleAdjointOperation(const Real& value, const IndexType lhsIndex, const Real* jacobies, const IndexType* rhsIndices, const int size);

    /**
     * @brief Pre definition of the the expression handles.
     *
     * @tparam AdjointData  The type of the adjoint data.
     * @tparam        Real  The type for the real values.
     * @tparam   IndexType  The index types for the management.
     */
    template<typename AdjointData, typename Real, typename IndexType> class ExpressionHandle;

    /**
     * @brief A function that is called for every statement that is written on the
     *        primal value tapes.
     *
     * The function can be used to extract information from the taping process.
     *
     * It can be set with the preprocessor macro CODI_AdjointHandle=<true/false>
     *
     * @param[in]          value  The primal value of the statement.
     * @param[in]       lhsIndex  The index on the left hand side of the statement, that
     *                            will be recorded.
     * @param[in]         handle  The handle that describes the whole expression that will be recorded on the tape.
     * @param[in] passiveActives  The number of active real variables that are passive in the statement (e.g. index == 0)
     * @param[in]      constants  The constant values that will be stored on the tape. The number
     *                            is the constant variable count from the expression and the passiveActives number.
     * @param[in]     rhsIndices  The pointer to the array of the indices that will be stored
     *                            for the statement. rhsIndices[0] is the first argument, rhsIndices[1]
     *                            the second, etc.
     * @param[in]      primalVec  The global vector of the primal variables.
     *
     * @tparam        Real  The type of the floating point values that are used in the tape.
     * @tparam PassiveReal  The type of the passive floating point values that are used in the tape.
     * @tparam   IndexType  The type of the indices that are used in the tape.
     */
    template<typename Real, typename PassiveReal, typename IndexType>
    void handleAdjointOperation(const Real& value, const IndexType lhsIndex, const ExpressionHandle<Real*, Real, IndexType>* handle, const StatementInt& passiveActives, const PassiveReal* constants, const IndexType* rhsIndices, const Real* primalVec);

    /**
     * @brief A function that is called for every statement that is evaluated in the forward tape.
     *
     * The function can be used to extract information from the taping process.
     *
     * It can be set with the preprocessor macro CODI_AdjointHandle=<true/false>
     *
     * @param[in]   value  The primal value of the statement.
     * @param[in] tangent  The tangent value of the statement.
     *
     * @tparam        Real  The type of the floating point values that are used in the tape.
     * @tparam TangentReal  The type of the tangent value that is used in the tape.
     */
    template<typename Real, typename TangentReal>
    void handleTangentOperation(const Real& value, const TangentReal& tangent);
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
       *
       * It can be set with the preprocessor macro CODI_EnableAssert=<true/false>
       */
      #define codiAssert(x) codi::checkAndOutputAssert(x, CODI_TO_STRING(x), __PRETTY_FUNCTION__, __FILE__, __LINE__)
    #else
      /**
       * @brief The assert function for CoDiPack it can be enabled with the preprocessor macro CODI_EnableAssert=true
       *
       * It can be set with the preprocessor macro CODI_EnableAssert=<true/false>
       *
       * @param x The expression that is checked in the assert.
       */
      #define codiAssert(x) /* disabled by CODI_EnableAssert */
    #endif
  #endif
}
