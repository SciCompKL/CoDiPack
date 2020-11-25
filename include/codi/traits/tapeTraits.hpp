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

  namespace TapeTraits {

    template<typename Tape>
    using IsForwardTape = std::is_base_of<ForwardEvaluation<typename Tape::Real, typename Tape::Gradient>, Tape>;

#if CODI_IS_CPP14
    template<typename Tape>
    bool constexpr isForwardTape = IsForwardTape<Tape>::value;
#endif

    template<typename Tape>
    using enableIfForwardTape = typename std::enable_if<IsForwardTape<Tape>::value>::type;

    template<typename Tape>
    using IsPrimalValueTape = std::is_base_of<PrimalValueBaseTape<typename Tape::TapeTypes, Tape>, Tape>;

#if CODI_IS_CPP14
    template<typename Tape>
    bool constexpr isPrimalValueTape = IsPrimalValueTape<Tape>::value;
#endif

    template<typename Tape>
    using enableIfPrimalValueTape = typename std::enable_if<IsPrimalValueTape<Tape>::value>::type;

    template<typename Tape>
    using IsJacobianTape = std::is_base_of<JacobianBaseTape<typename Tape::TapeTypes, Tape>, Tape>;

#if CODI_IS_CPP14
    template<typename Tape>
    bool constexpr isJacobianTape = IsJacobianTape<Tape>::value;
#endif

    template<typename Tape>
    using enableIfJacobianTape = typename std::enable_if<IsJacobianTape<Tape>::value>::type;

    template<typename Tape>
    using IsReverseTape = std::integral_constant<bool,
           IsJacobianTape<Tape>::value
        || IsPrimalValueTape<Tape>::value
      >;

#if CODI_IS_CPP14
    template<typename Tape>
    bool constexpr isReverseTape = IsReverseTape<Tape>::value;
#endif

    template<typename Tape>
    using enableIfReverseTape = typename std::enable_if<IsReverseTape<Tape>::value>::type;
  }
}
