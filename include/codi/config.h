#pragma once

#include <stdint.h>
#include <stddef.h>

/** @file */

/** \copydoc codi::Namespace */
namespace codi {

  // Can be enabled to allow for auto completion in IDEs.
  #ifndef CODI_IDE
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

  namespace Config {

    /*******************************************************************************/
    /// @name Type and compile time value declarations
    /// @{

    #ifndef CODI_ChunkSize
      #define CODI_ChunkSize 2097152
    #endif
    /// Default size of chunks (ChunkBase) used in ChunkedData in reverse tape implementations.
    static size_t ChunkSize = CODI_ChunkSize;
    #undef CODI_ChunkSize

    /// Type for the number of arguments in statements
    using ArgumentSize = uint8_t;

    /// Maximum number of arguments in a statement
    size_t constexpr MaxArgumentSize = 255;

    /// Tag for statements that are inputs. Used in linear index management context.
    size_t constexpr StatementInputTag = 255;

    #ifndef CODI_SmallChunkSize
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
      #define CODI_EnableMPI false
    #endif
    /// Add MPI and MeDiPack specific headers
    bool constexpr EnableMPI = CODI_EnableMPI;
    // Do not undefine

    #ifndef CODI_CheckExpressionArguments
      #define CODI_CheckExpressionArguments false
    #endif
    /// Check for invalid arguments to expressions like division by zero.
    bool constexpr CheckExpressionArguments = CODI_CheckExpressionArguments;
    #undef CODI_CheckExpressionArguments

    #ifndef CODI_CheckJacobiIsZero
      #define CODI_CheckJacobiIsZero true
    #endif
    /// Ignore Jacobian that are zero in Jacobian based tapes.
    bool constexpr CheckJacobiIsZero = CODI_CheckJacobiIsZero;
    #undef CODI_CheckJacobiIsZero

    #ifndef CODI_CheckTapeActivity
      #define CODI_CheckTapeActivity true
    #endif
    /// Makes it possible to ignore certain code parts. If turned of everything will be recorded.
    bool constexpr CheckTapeActivity = CODI_CheckTapeActivity;
    #undef CODI_CheckTapeActivity

    #ifndef CODI_CheckZeroIndex
      #define CODI_CheckZeroIndex true
    #endif
    /// Ignore active types that are not dependent on any input value in Jacobian tapes.
    bool constexpr CheckZeroIndex = CODI_CheckZeroIndex;
    #undef CODI_CheckZeroIndex

    #ifndef CODI_CopyOptimization
      #define CODI_CopyOptimization true
    #endif
    /// Do not store copy statements like a = b; if the identity handler allows it.
    bool constexpr CopyOptimization = CODI_CopyOptimization;
    #undef CODI_CopyOptimization

    #ifndef CODI_RemoveDuplicateJacobianArguments
      #define CODI_RemoveDuplicateJacobianArguments 0
    #endif
    /// Extra pass in Jacobian tapes that combines arguments with the same identifier.
    bool constexpr RemoveDuplicateJacobianArguments = CODI_RemoveDuplicateJacobianArguments;
    // Do not undefine

    #ifndef CODI_IgnoreInvalidJacobies
      #define CODI_IgnoreInvalidJacobies false
    #endif
    /// Ignore invalid Jacobians like NaN or Inf
    bool constexpr IgnoreInvalidJacobies = CODI_IgnoreInvalidJacobies;
    #undef CODI_IgnoreInvalidJacobies

    #ifndef CODI_OverflowCheck
      #define CODI_OverflowCheck true
    #endif
    /// Check in the index manager if an overflow occurred.
    bool constexpr OverflowCheck = CODI_OverflowCheck;
    #undef CODI_OverflowCheck

    #ifndef CODI_SkipZeroAdjointEvaluation
      #define CODI_SkipZeroAdjointEvaluation true
    #endif
    /// Do not perform a reverse evaluation of a statement if the seeding adjoint is zero.
    bool constexpr SkipZeroAdjointEvaluation = CODI_SkipZeroAdjointEvaluation;
    #undef CODI_SkipZeroAdjointEvaluation

    #ifndef CODI_SortIndicesOnReset
      #define CODI_SortIndicesOnReset true
    #endif
    /// Reuse index tapes will sort their indices on a reset.
    bool constexpr SortIndicesOnReset = CODI_SortIndicesOnReset;
    #undef CODI_SortIndicesOnReset

    #ifndef CODI_VariableAdjointInterfaceInPrimalTapes
      #define CODI_VariableAdjointInterfaceInPrimalTapes 0
    #endif
    /// Allow custom adjoint vector in primal values tapes.
    bool constexpr VariableAdjointInterfaceInPrimalTapes = CODI_VariableAdjointInterfaceInPrimalTapes;
    #if CODI_VariableAdjointInterfaceInPrimalTapes
      #define ADJOINT_VECTOR_TYPE VectorAccessInterface<Real, Identifier>
    #else
      #define ADJOINT_VECTOR_TYPE Gradient
    #endif
    // Do not undefine

    /// @}
    /*******************************************************************************/
    /// @name Macro definitions
    /// @{

    #ifndef CODI_AvoidedInlines
      #define CODI_AvoidedInlines 1
    #endif
    #if CODI_AvoidedInlines
      #if defined(_MSC_VER)
        #define CODI_NO_INLINE __declspec(noinline)
      #else
        #define CODI_NO_INLINE __attribute__((noinline))
      #endif
    #else
      #define CODI_NO_INLINE /* no avoiding of inline defined */
    #endif
    /// Do not inline functions like evaluate()
    bool constexpr AvoidedInlines = CODI_AvoidedInlines;
    #undef CODI_AvoidedInlines

    #ifndef CODI_EnableAssert
      #define CODI_EnableAssert false
    #endif
    #ifndef codiAssert // Can also be defined by the user before including codi.hpp
      #if CODI_EnableAssert
        #define codiAssert(x) ::codi::checkAndOutputAssert(x, CODI_TO_STRING(x), __PRETTY_FUNCTION__, __FILE__, __LINE__)
      #else
        #define codiAssert(x) /* disabled by CODI_EnableAssert */
      #endif
    #endif
    /// Enables asserts in CoDiPack for consistency checking.
    bool constexpr EnableAssert = CODI_EnableAssert;
    #undef CODI_EnableAssert

    #ifndef CODI_ForcedInlines
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
      #define CODI_INLINE inline
    #endif
    /// Force inlining instead of using the heuristics from the compiler.
    bool constexpr ForcedInlines = CODI_ForcedInlines;
    #undef CODI_ForcedInlines

    /// @}
  }
}
