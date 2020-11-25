#pragma once

#include <type_traits>

#include "../../aux/macros.hpp"
#include "../../config.h"

/** \copydoc codi::Namespace */
namespace codi {

  /// Enable if abbreviation for std::is_base_of
  template<typename Base, typename Impl>
  using enable_if_base_of = std::enable_if<std::is_base_of<Base, Impl>::value>;

  /// Enable if abbreviation for std::is_same
  template<typename T1, typename T2>
  using enable_if_same = std::enable_if<std::is_same<T1, T2>::value>;
}
