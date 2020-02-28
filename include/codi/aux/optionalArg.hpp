#pragma once

#include "../config.h"

/** \copydoc codi::Namespace */
namespace codi {

  template<typename T>
  struct OptionalArg {
      static T value;
  };

  template<typename T>
  T OptionalArg<T>::value = {};
}
