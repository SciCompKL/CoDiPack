/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2023 Chair for Scientific Computing (SciComp), University of Kaiserslautern-Landau
 * Homepage: http://www.scicomp.uni-kl.de
 * Contact:  Prof. Nicolas R. Gauger (codi@scicomp.uni-kl.de)
 *
 * Lead developers: Max Sagebaum, Johannes Blühdorn (SciComp, University of Kaiserslautern-Landau)
 *
 * This file is part of CoDiPack (http://www.scicomp.uni-kl.de/software/codi).
 *
 * CoDiPack is free software: you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * CoDiPack is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty
 * of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 * See the GNU General Public License for more details.
 * You should have received a copy of the GNU
 * General Public License along with CoDiPack.
 * If not, see <http://www.gnu.org/licenses/>.
 *
 * For other licensing options please contact us.
 *
 * Authors:
 *  - SciComp, University of Kaiserslautern-Landau:
 *    - Max Sagebaum
 *    - Johannes Blühdorn
 *    - Former members:
 *      - Tim Albring
 */
#pragma once

#include <type_traits>

#include "../config.h"
#include "../misc/macros.hpp"
#include "../tapes/jacobianBaseTape.hpp"
#include "../tapes/primalValueBaseTape.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  template<typename T_Real, typename T_Gradient>
  struct ForwardEvaluation;

  /// Traits for everything that can be a CoDiPack tape usually the template argument of codi::ActiveType.
  /// Possible types are codi::JacobianLinearTape, codi::ForwardEvaluation, codi::PrimalValueReuseTape, etc..
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
    using IsReverseTape = std::integral_constant<bool, IsJacobianTape<Tape>::value || IsPrimalValueTape<Tape>::value>;

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
