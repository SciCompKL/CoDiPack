#pragma once

#include "../config.h"

/** \copydoc codi::Namespace */
namespace codi {


  /*******************************************************************************
   * Section: Default type declarations
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
  #define ANY void

  /*
   * Used in interface declarations for types that have to be defined in the specializations.
   */
  #define UNDEFINED void

  /*
   * Expand template types in preprocessor macros.
   */
  #define TEMPLATE(...) __VA_ARGS__
}
