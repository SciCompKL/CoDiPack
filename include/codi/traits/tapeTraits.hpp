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

    /*******************************************************************************/
    /// @name Detection of specific real value types
    /// @{

    /// If the tape inherits from ForwardEvaluation.
    template<typename Tape>
    using IsForwardTape = std::is_base_of<ForwardEvaluation<typename Tape::Real, typename Tape::Gradient>, Tape>;

#if CODI_IS_CPP14
    /// Value entry of IsForwardTape
    template<typename Tape>
    bool constexpr isForwardTape = IsForwardTape<Tape>::value;
#endif

    /// Enable if wrapper for IsForwardTape
    template<typename Tape>
    using EnableIfForwardTape = typename std::enable_if<IsForwardTape<Tape>::value>::type;

    /// If the tape inherits from PrimalValueBaseTape.
    template<typename Tape>
    using IsPrimalValueTape = std::is_base_of<PrimalValueBaseTape<typename Tape::TapeTypes, Tape>, Tape>;

#if CODI_IS_CPP14
    /// Value entry of IsPrimalValueTape
    template<typename Tape>
    bool constexpr isPrimalValueTape = IsPrimalValueTape<Tape>::value;
#endif

    /// Enable if wrapper for IsPrimalValueTape
    template<typename Tape>
    using EnableIfPrimalValueTape = typename std::enable_if<IsPrimalValueTape<Tape>::value>::type;

    /// If the tape inherits from JacobianBaseTape.
    template<typename Tape>
    using IsJacobianTape = std::is_base_of<JacobianBaseTape<typename Tape::TapeTypes, Tape>, Tape>;

#if CODI_IS_CPP14
    /// Value entry of IsJacobianTape
    template<typename Tape>
    bool constexpr isJacobianTape = IsJacobianTape<Tape>::value;
#endif

    /// Enable if wrapper for IsJacobianTape
    template<typename Tape>
    using EnableIfJacobianTape = typename std::enable_if<IsJacobianTape<Tape>::value>::type;

    /// If the tape inherits either from JacobianBaseTape or PrimalValueBaseTape.
    template<typename Tape>
    using IsReverseTape = std::integral_constant<bool,
           IsJacobianTape<Tape>::value
        || IsPrimalValueTape<Tape>::value
      >;

#if CODI_IS_CPP14
    /// Value entry of IsReverseTape
    template<typename Tape>
    bool constexpr isReverseTape = IsReverseTape<Tape>::value;
#endif

    /// Enable if wrapper for IsReverseTape
    template<typename Tape>
    using EnableIfReverseTape = typename std::enable_if<IsReverseTape<Tape>::value>::type;

    /// @}
  }
}
