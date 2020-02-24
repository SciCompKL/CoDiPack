#pragma once

#include "../config.h"

/** \copydoc codi::Namespace */
namespace codi {

  /*
   * Expand template types in preprocessor macros.
   */
  #define TEMPLATE(...) __VA_ARGS__

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
}
