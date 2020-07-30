#pragma once

#include "../config.h"

/** \copydoc codi::Namespace */
namespace codi {

  template<typename ...Args>
  void CODI_UNUSED(Args const& ... ) {}

  /* used in a constant context, where using CODI_UNUSED spoils the constantness */
  #define CODI_UNUSED_ARG(arg) /* arg */

  #define CODI_ENABLE_CHECK(option, condition) if(!(option) || (condition))

  #define CODI_TO_STRING2(expression) #expression

  #define CODI_TO_STRING(expression) CODI_TO_STRING2(expression)

  /*******************************************************************************
   * Section: Default template type declarations
   *
   * Description: TODO
   *
   */

  /*
   * CODI_IDE can be define to use the default declaration of typenames. This enables autocompletion in the IDEs.
   *
   * Every using declartion in all CoDiPack classes should declare its variables as:
   *  using TYPE = CODI_DECLARE_DEFAULT(_TYPE, Default);
   */
  #if CODI_IDE
    #define CODI_DECLARE_DEFAULT(Type, Default) Default
  #else
    #define CODI_DECLARE_DEFAULT(Type, Default) Type
  #endif

  /*
   * Used in default declarations of expression templates.
   */
  #define CODI_ANY int

  /*
   * Expand template types in preprocessor macros.
   */
  #define CODI_TEMPLATE(...) __VA_ARGS__

  /*
   * Used in interface declarations for types that have to be defined in the specializations.
   */
  #define CODI_UNDEFINED void

  /*
   * Used in interface declarations for variables that have to be defined in the specializations.
   */
  #define CODI_UNDEFINED_VALUE false

  template<typename First, typename... Tail>
  struct CODI_UNION : public First, public CODI_UNION<Tail...> {};

  template<typename First>
  struct CODI_UNION<First> : public First {};

  #define CODI_WRAP_FUNCTION(NAME, FUNC) \
    struct NAME { \
      template<typename ... Args> \
      void operator()(Args&& ... args) const { \
        FUNC(std::forward<Args>(args)...); \
      } \
    }

  #define CODI_WRAP_FUNCTION_TEMPLATE(NAME, FUNC) \
    template<typename ... TT> \
    struct NAME { \
      template<typename ... Args> \
      void operator()(Args&& ... args) const { \
        FUNC<TT...>(std::forward<Args>(args)...); \
      } \
    }
}
