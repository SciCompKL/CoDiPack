#pragma once

/** \copydoc codi::Namespace */
namespace codi {

  #ifndef IDE
    #define IDE 0
  #endif

  /**
   * @brief TODO: CoDiPack namespace documentation.
   */
  struct Namespace {};

  struct Config {
    public:

      #ifndef CODI_CheckExpressionArguments
        #define CODI_CheckExpressionArguments false
      #endif
      static bool constexpr CheckExpressionArguments = CODI_CheckExpressionArguments;
      #undef CODI_CheckExpressionArguments

      #ifndef CODI_IgnoreInvalidJacobies
        #define CODI_IgnoreInvalidJacobies false
      #endif
      static bool constexpr IgnoreInvalidJacobies = CODI_IgnoreInvalidJacobies;
      #undef CODI_OptIgnoreInvalidJacobies


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
      static bool constexpr IsForcedInlines = CODI_UseForcedInlines;
      #undef CODI_UseForcedInlines

  };
}
