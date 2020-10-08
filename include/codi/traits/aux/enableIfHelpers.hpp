#pragma once

#include <type_traits>

#include "../../aux/macros.hpp"
#include "../../config.h"

/** \copydoc codi::Namespace */
namespace codi {

  template<typename Base, typename Impl>
  using enable_if_base_of = std::enable_if<std::is_base_of<Base, Impl>::value>;

  template<typename T1, typename T2>
  using enable_if_same = std::enable_if<std::is_same<T1, T2>::value>;
}
