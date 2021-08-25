#pragma once

#include <stddef.h>
#include <stdint.h>

#include "aux/exceptions.hpp"

/** @file */

/** \copydoc codi::Namespace */
namespace codi {

#ifndef CODI_IDE
  /// Can be enabled to allow for auto completion in IDEs.
  #define CODI_IDE 0
#endif

  /**
   * @brief CoDiPack - Code Differentiation Package
   *
   * Web: https://www.scicomp.uni-kl.de/software/codi
   * Git: https://github.com/scicompkl/codipack
   * Doc: https://www.scicomp.uni-kl.de/codi
   */
  struct Namespace {};

  /// Configuration options for CoDiPack.
  namespace Config {

    /*******************************************************************************/
    /// @name Type and compile time value declarations
    /// @{

#ifndef CODI_ChunkSize
  /// See codi::Config::ChunkSize
  #define CODI_ChunkSize 2097152
#endif
    /// Default size of chunks (ChunkBase) used in ChunkedData in reverse tape implementations.
    size_t constexpr ChunkSize = CODI_ChunkSize;
#undef CODI_ChunkSize

    /// Type for the number of arguments in statements
    using ArgumentSize = uint8_t;

    /// Maximum number of arguments in a statement
    size_t constexpr MaxArgumentSize = 254;

    /// Tag for statements that are inputs. Used in linear index management context.
    size_t constexpr StatementInputTag = 255;

#ifndef CODI_SmallChunkSize
  /// See codi::Config::SmallChunkSize
  #define CODI_SmallChunkSize 32768
#endif
    /// Default smaller size of chunks (ChunkBase) used in ChunkedData in reverse tape implementations.
    size_t constexpr SmallChunkSize = CODI_SmallChunkSize;
#undef CODI_SmallChunkSize

    /// @}
    /*******************************************************************************/
    /// @name Compile time flags
    /// @{

#ifndef CODI_EnableMPI
  /// See codi::Config::EnableMPI
  #define CODI_EnableMPI false
#endif
    /// Add MPI and MeDiPack specific headers
    bool constexpr EnableMPI = CODI_EnableMPI;
    // Do not undefine.

#ifndef CODI_CheckExpressionArguments
  /// See codi::Config::CheckExpressionArguments
  #define CODI_CheckExpressionArguments false
#endif
    /// Check for invalid arguments to expressions like division by zero.
    bool constexpr CheckExpressionArguments = CODI_CheckExpressionArguments;
#undef CODI_CheckExpressionArguments

#ifndef CODI_CheckJacobianIsZero
  /// See codi::Config::CheckJacobianIsZero
  #define CODI_CheckJacobianIsZero true
#endif
    /// Ignore Jacobians that are zero in Jacobian based tapes.
    bool constexpr CheckJacobianIsZero = CODI_CheckJacobianIsZero;
#undef CODI_CheckJacobianIsZero

#ifndef CODI_CheckTapeActivity
  /// See codi::Config::CheckTapeActivity
  #define CODI_CheckTapeActivity true
#endif
    /// Makes it possible to ignore certain code parts. If turned of everything will be recorded.
    bool constexpr CheckTapeActivity = CODI_CheckTapeActivity;
#undef CODI_CheckTapeActivity

#ifndef CODI_CheckZeroIndex
  /// See codi::Config::CheckZeroIndex
  #define CODI_CheckZeroIndex true
#endif
    /// Ignore active types that are not dependent on any input value in Jacobian tapes.
    bool constexpr CheckZeroIndex = CODI_CheckZeroIndex;
#undef CODI_CheckZeroIndex

#ifndef CODI_CopyOptimization
  /// See codi::Config::CopyOptimization
  #define CODI_CopyOptimization true
#endif
    /// Do not store copy statements like a = b; if the identity handler allows it.
    bool constexpr CopyOptimization = CODI_CopyOptimization;
#undef CODI_CopyOptimization

#ifndef CODI_RemoveDuplicateJacobianArguments
  /// See codi::Config::RemoveDuplicateJacobianArguments
  #define CODI_RemoveDuplicateJacobianArguments 0
#endif
    /// Extra pass in Jacobian tapes that combines arguments with the same identifier.
    bool constexpr RemoveDuplicateJacobianArguments = CODI_RemoveDuplicateJacobianArguments;
    // Do not undefine.

#ifndef CODI_IgnoreInvalidJacobians
  /// See codi::Config::IgnoreInvalidJacobians
  #define CODI_IgnoreInvalidJacobians false
#endif
    /// Ignore invalid Jacobians like NaN or Inf
    bool constexpr IgnoreInvalidJacobians = CODI_IgnoreInvalidJacobians;
#undef CODI_IgnoreInvalidJacobians

#ifndef CODI_OverflowCheck
  /// See codi::Config::OverflowCheck
  #define CODI_OverflowCheck true
#endif
    /// Check in the index manager if an overflow occurred.
    bool constexpr OverflowCheck = CODI_OverflowCheck;
#undef CODI_OverflowCheck

#ifndef CODI_SkipZeroAdjointEvaluation
  /// See codi::Config::SkipZeroAdjointEvaluation
  #define CODI_SkipZeroAdjointEvaluation true
#endif
    /// Do not perform a reverse evaluation of a statement if the seeding adjoint is zero.
    bool constexpr SkipZeroAdjointEvaluation = CODI_SkipZeroAdjointEvaluation;
#undef CODI_SkipZeroAdjointEvaluation

#ifndef CODI_SortIndicesOnReset
  /// See codi::Config::SortIndicesOnReset
  #define CODI_SortIndicesOnReset true
#endif
    /// Reuse index tapes will sort their indices on a reset.
    bool constexpr SortIndicesOnReset = CODI_SortIndicesOnReset;
#undef CODI_SortIndicesOnReset

#ifndef CODI_VariableAdjointInterfaceInPrimalTapes
  /// See codi::Config::VariableAdjointInterfaceInPrimalTapes
  #define CODI_VariableAdjointInterfaceInPrimalTapes 0
#endif
    /// Allow custom adjoint vector in primal values tapes.
    bool constexpr VariableAdjointInterfaceInPrimalTapes = CODI_VariableAdjointInterfaceInPrimalTapes;
#if CODI_VariableAdjointInterfaceInPrimalTapes
  #define ADJOINT_VECTOR_TYPE VectorAccessInterface<Real, Identifier>
#else
  /// See codi::Config::VariableAdjointInterfaceInPrimalTapes
  #define ADJOINT_VECTOR_TYPE Gradient
#endif
    // Do not undefine.

#ifndef CODI_ReversalZeroesAdjoints
  /// See codi::Config::ReversalZeroesAdjoints
  #define CODI_ReversalZeroesAdjoints true
#endif
#if CODI_VariableAdjointInterfaceInPrimalTapes && !CODI_ReversalZeroesAdjoints
  #warning CODI_ReversalZeroesAdjoints == false is incompatible with CODI_VariableAdjointInterfaceInPrimalTapes == true.
#endif
    /// With a linear index management, control if adjoints are set to zero during reversal.
    bool constexpr ReversalZeroesAdjoints = CODI_ReversalZeroesAdjoints;
#undef CODI_ReversalZeroesAdjoints

    /// @}
    /*******************************************************************************/
    /// @name Macro definitions
    /// @{

#ifndef CODI_AvoidedInlines
  /// See codi::Config::AvoidedInlines
  #define CODI_AvoidedInlines 1
#endif
#if CODI_AvoidedInlines
  #if defined(_MSC_VER)
    #define CODI_NO_INLINE __declspec(noinline)
  #else
    #define CODI_NO_INLINE __attribute__((noinline))
  #endif
#else
  /// See codi::Config::AvoidedInlines
  #define CODI_NO_INLINE /* no avoiding of inline defined */
#endif
    /// Do not inline functions like evaluate()
    bool constexpr AvoidedInlines = CODI_AvoidedInlines;
#undef CODI_AvoidedInlines

#ifndef CODI_EnableAssert
  /// See codi::Config::EnableAssert
  #define CODI_EnableAssert false
#endif
#ifndef codiAssert  // Can also be defined by the user before including codi.hpp.
  #if CODI_EnableAssert
    #define codiAssert(x) ::codi::checkAndOutputAssert(x, CODI_TO_STRING(x), __PRETTY_FUNCTION__, __FILE__, __LINE__)
  #else
    /// See codi::Config::EnableAssert
    #define codiAssert(x) /* disabled by CODI_EnableAssert */
  #endif
#endif
    /// Enables asserts in CoDiPack for consistency checking.
    bool constexpr EnableAssert = CODI_EnableAssert;
#undef CODI_EnableAssert

#ifndef CODI_ForcedInlines
  /// See codi::Config::ForcedInlines
  #define CODI_ForcedInlines 0
#endif
#if CODI_ForcedInlines
  #if defined(__INTEL_COMPILER) | defined(_MSC_VER)
    #define CODI_INLINE __forceinline
  #elif defined(__GNUC__)
    #define CODI_INLINE inline __attribute__((always_inline))
  #else
    #warning Could not determine compiler for forced inline definitions. Using inline.
    #define CODI_INLINE inline
  #endif
#else
  /// See codi::Config::ForcedInlines
  #define CODI_INLINE inline
#endif
    /// Force inlining instead of using the heuristics from the compiler.
    bool constexpr ForcedInlines = CODI_ForcedInlines;
#undef CODI_ForcedInlines

    /// @}
  }
}
