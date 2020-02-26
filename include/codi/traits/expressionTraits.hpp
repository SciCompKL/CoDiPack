#pragma once

#include <type_traits>

#include "../aux/macros.h"
#include "../config.h"

/** \copydoc codi::Namespace */
namespace codi {

  template<typename _Real, typename _Gradient, typename _Tape, typename _Impl>
  struct LhsExpressionInterface;

  template<typename Real, typename Gradient, typename Tape, typename Impl>
  using isLhsExpression = std::is_base_of<LhsExpressionInterface<Real, Gradient, Tape, Impl>, Impl>;

  template<typename Real, typename Gradient, typename Tape, typename Impl>
  using enableIfLhsExpression = typename std::enable_if<isLhsExpression<Real, Gradient, Tape, Impl>::value>::type;
}
