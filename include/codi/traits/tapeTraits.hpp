#pragma once

#include <type_traits>

#include "../aux/macros.hpp"
#include "../config.h"

#include "../tapes/jacobianBaseTape.hpp"
#include "../tapes/primalValueBaseTape.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  template<typename _Real, typename _Gradient>
  struct ForwardEvaluation;

  template<typename Tape>
  using isForwardTape = std::is_base_of<ForwardEvaluation<typename Tape::Real, typename Tape::Gradient>, Tape>;

  template<typename Tape>
  using enableIfForwardTape = typename std::enable_if<isForwardTape<Tape>::value>::type;

  template<typename Tape>
  using isPrimalValueTape = std::is_base_of<PrimalValueBaseTape<typename Tape::TapeTypes, Tape>, Tape>;

  template<typename Tape>
  using enableIfPrimalValueTape = typename std::enable_if<isPrimalValueTape<Tape>::value>::type;

  template<typename Tape>
  using isJacobianTape = std::is_base_of<JacobianBaseTape<typename Tape::TapeTypes, Tape>, Tape>;

  template<typename Tape>
  using enableIfJacobianTape = typename std::enable_if<isJacobianTape<Tape>::value>::type;

  template<typename Tape>
  using isReverseTape = std::integral_constant<bool,
         isJacobianTape<Tape>::value
      || isPrimalValueTape<Tape>::value
    >;

  template<typename Tape>
  using enableIfReverseTape = typename std::enable_if<isReverseTape<Tape>::value>::type;
}
