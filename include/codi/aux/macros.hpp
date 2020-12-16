#pragma once

#include "../config.h"

/** \copydoc codi::Namespace */
namespace codi {

  /// Disable unused warnings for arbitrary number of arguments.
  template<typename ...Args>
  void CODI_UNUSED(Args const& ... ) {}

  /// used in a constant context, where using CODI_UNUSED spoils the constantness */
  #define CODI_UNUSED_ARG(arg) /* arg */

  /// Enable the if only if the option is true otherwise the code block is always evaluated
  #define CODI_ENABLE_CHECK(option, condition) if(!(option) || (condition))

  /// Conversion macro
  #define CODI_TO_STRING2(expression) #expression

  /// Conversion macro
  #define CODI_TO_STRING(expression) CODI_TO_STRING2(expression)

  /// Check for CPP 14 standard
  #define CODI_IS_CPP14 (201402L <= __cplusplus)

  /// Check for CPP 17 standard
  #define CODI_IS_CPP17 (201703L <= __cplusplus)


  /*******************************************************************************/
  /** @name Default template type declarations
   *
   * Description: These template are used to employ the design guide line for the default definitions of template
   *              arguments. See \ref TemplateDeclaration for Details.
   * @{
   */

  /**
   * CODI_IDE can be defined to use the default declaration of type names. This enables auto completion in the IDEs.
   *
   * Every using declaration in all CoDiPack classes should declare its variables as:
   *  using TYPE = CODI_DECLARE_DEFAULT(_TYPE, Default);
   */
  #if CODI_IDE
    #define CODI_DECLARE_DEFAULT(Type, Default) Default
  #else
    #define CODI_DECLARE_DEFAULT(Type, Default) Type
  #endif

  /// Abbreviation for CODI_DECLARE_DEFAULT
  #define CODI_DD(Type, Default) CODI_DECLARE_DEFAULT(Type, Default)

  /// Used in default declarations of expression templates.
  #define CODI_ANY int

  /// Used in interface declarations to indicate the type of the implementing class.
  #define CODI_IMPLEMENTATION int

  /// Expand template types in preprocessor macros.
  #define CODI_TEMPLATE(...) __VA_ARGS__

  /// Abbreviation for CODI_TEMPLATE
  #define CODI_T(...) CODI_TEMPLATE(__VA_ARGS__)

  /// Used in interface declarations for types that have to be defined in the specializations.
  #define CODI_UNDEFINED void

  /// Used in interface declarations for variables that have to be defined in the specializations.
  #define CODI_UNDEFINED_VALUE false

  /// Creates a union of interface definitions
  template<typename First, typename... Tail>
  struct CODI_UNION : public First, public CODI_UNION<Tail...> {};

  template<typename First>
  struct CODI_UNION<First> : public First {};

  /// @}

  /// Wrap a function in a function object. Used for speed optimizations.
  #define CODI_WRAP_FUNCTION(NAME, FUNC) \
    struct NAME { \
      /** Empty */ \
      template<typename ... Args> \
      void operator()(Args&& ... args) const { \
        FUNC(std::forward<Args>(args)...); \
      } \
  }

  /// Wrap a function in a function object. Used for speed optimizations.
  #define CODI_WRAP_FUNCTION_TEMPLATE(NAME, FUNC) \
    template<typename ... TT> \
    struct NAME { \
      /** Empty */ \
      template<typename ... Args> \
      void operator()(Args&& ... args) const { \
        FUNC<TT...>(std::forward<Args>(args)...); \
      } \
    }
}
