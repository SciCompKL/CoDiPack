#pragma once

#include "../config.h"

/** \copydoc codi::Namespace */
namespace codi {

  template<typename ...Args>
  void CODI_UNUSED(Args const& ... ) {}

  #define CODI_UNUSED_ARG(arg) /* arg */

  #define ENABLE_CHECK(option, condition) if(!(option) || (condition))

  #define CODI_TO_STRING2(expression) #expression

  #define CODI_TO_STRING(expression) CODI_TO_STRING2(expression)

  /*******************************************************************************
   * Section: Default template type declarations
   *
   * Description: TODO
   *
   */

  /*
   * IDE can be define to use the default declaration of typenames. This enables autocompletion in the IDEs.
   *
   * Every using declartion in all CoDiPack classes should declare its variables as:
   *  using TYPE = DECLARE_DEFAULT(_TYPE, Default);
   */
  #if IDE
    #define DECLARE_DEFAULT(Type, Default) Default
  #else
    #define DECLARE_DEFAULT(Type, Default) Type
  #endif

  /*
   * Used in default declarations of expression templates.
   */
  #define ANY int

  /*
   * Used in interface declarations for types that have to be defined in the specializations.
   */
  #define CODI_UNDEFINED void

  /*
   * Used in interface declarations for variables that have to be defined in the specializations.
   */
  #define CODI_UNDEFINED_VALUE false

  /*
   * Expand template types in preprocessor macros.
   */
  #define TEMPLATE(...) __VA_ARGS__

  #define WRAP_FUNCTION(NAME, FUNC) \
    struct NAME { \
      template<typename ... Args> \
      void operator()(Args&& ... args) const { \
        FUNC(std::forward<Args>(args)...); \
      } \
    }

#define WRAP_FUNCTION_TEMPLATE(NAME, FUNC) \
  template<typename ... TT> \
  struct NAME { \
    template<typename ... Args> \
    void operator()(Args&& ... args) const { \
      FUNC<TT...>(std::forward<Args>(args)...); \
    } \
  }
}
