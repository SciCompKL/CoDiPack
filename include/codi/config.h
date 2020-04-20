#pragma once

#include <stdint.h>
#include <stddef.h>

/** \copydoc codi::Namespace */
namespace codi {

  #ifndef IDE
    #define IDE 0
  #endif

  /**
   * @brief TODO: CoDiPack namespace documentation.
   */
  struct Namespace {};

  namespace Config {

      /*******************************************************************************
       * Section: Type and compile time value declarations
       *
       * Description: TODO
       *
       */

      #ifndef CODI_ChunkSize
        #define CODI_ChunkSize 2097152
      #endif
      static size_t ChunkSize = CODI_ChunkSize;
      #undef CODI_ChunkSize

      using ArgumentSize = uint8_t;
      size_t constexpr MaxArgumentSize = 254;

      size_t constexpr StatementInputTag = 255;

      #ifndef CODI_SmallChunkSize
        #define CODI_SmallChunkSize 32768
      #endif
      size_t constexpr SmallChunkSize = CODI_SmallChunkSize;
      #undef CODI_SmallChunkSize

      /*******************************************************************************
       * Section: Compile time flags
       *
       * Description: TODO
       *
       */

      #ifndef CODI_AssignOptimization
        #define CODI_AssignOptimization true
      #endif
      const bool AssignOptimization = CODI_AssignOptimization;
      #undef CODI_AssignOptimization

      #ifndef CODI_CheckExpressionArguments
        #define CODI_CheckExpressionArguments false
      #endif
      bool constexpr CheckExpressionArguments = CODI_CheckExpressionArguments;
      #undef CODI_CheckExpressionArguments

      #ifndef CODI_CheckJacobiIsZero
        #define CODI_CheckJacobiIsZero true
      #endif
      const bool CheckJacobiIsZero = CODI_CheckJacobiIsZero;
      #undef CODI_CheckJacobiIsZero

      #ifndef CODI_CheckTapeActivity
        #define CODI_CheckTapeActivity true
      #endif
      const bool CheckTapeActivity = CODI_CheckTapeActivity;
      #undef CODI_CheckTapeActivity

      #ifndef CODI_CheckZeroIndex
        #define CODI_CheckZeroIndex true
      #endif
      const bool CheckZeroIndex = CODI_CheckZeroIndex;
      #undef CODI_CheckZeroIndex

      #ifndef CODI_CombineJacobianArguments
        #define CODI_CombineJacobianArguments 0
      #endif
      bool constexpr CombineJacobianArguments = CODI_CombineJacobianArguments;

      #ifndef CODI_IgnoreInvalidJacobies
        #define CODI_IgnoreInvalidJacobies false
      #endif
      bool constexpr IgnoreInvalidJacobies = CODI_IgnoreInvalidJacobies;
      #undef CODI_IgnoreInvalidJacobies

      #ifndef CODI_OverflowCheck
        #define CODI_OverflowCheck true
      #endif
      const bool OverflowCheck = CODI_OverflowCheck;
      #undef CODI_OverflowCheck

      #ifndef CODI_SkipZeroAdjointEvaluation
        #define CODI_SkipZeroAdjointEvaluation true
      #endif
      const bool SkipZeroAdjointEvaluation = CODI_SkipZeroAdjointEvaluation;
      #undef CODI_SkipZeroAdjointEvaluation

      #ifndef CODI_SortIndicesOnReset
        #define CODI_SortIndicesOnReset true
      #endif
      const bool SortIndicesOnReset = CODI_SortIndicesOnReset;
      #undef CODI_SortIndicesOnReset

      #ifndef CODI_VariableAdjointInterfaceInPrimalTapes
        #define CODI_VariableAdjointInterfaceInPrimalTapes 0
      #endif
      #if CODI_VariableAdjointInterfaceInPrimalTapes
        #define ADJOINT_VECTOR_TYPE VectorAccessInterface<Real, Identifier>
      #else
        #define ADJOINT_VECTOR_TYPE Gradient
      #endif

      /*******************************************************************************
       * Section: Macro definitions
       *
       * Description: TODO
       *
       */

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
      bool constexpr ForcedInlines = CODI_ForcedInlines;
      #undef CODI_ForcedInlines
  }
}
