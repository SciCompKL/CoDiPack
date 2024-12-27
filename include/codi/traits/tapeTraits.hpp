/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2025 Chair for Scientific Computing (SciComp), University of Kaiserslautern-Landau
 * Homepage: http://scicomp.rptu.de
 * Contact:  Prof. Nicolas R. Gauger (codi@scicomp.uni-kl.de)
 *
 * Lead developers: Max Sagebaum, Johannes Blühdorn (SciComp, University of Kaiserslautern-Landau)
 *
 * This file is part of CoDiPack (http://scicomp.rptu.de/software/codi).
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
#include "../tapes/interfaces/editingTapeInterface.hpp"
#include "misc/enableIfHelpers.hpp"

/** \copydoc codi::Namespace */
namespace codi {
  template<typename Tape, typename Impl>
  struct PrimalValueBaseTape;

  template<typename Tape, typename Impl>
  struct JacobianBaseTape;

  template<typename T_Real, typename T_Gradient>
  struct ForwardEvaluation;

  template<typename T_Real, typename T_Tag>
  struct TagTapeForward;

  template<typename T_Real, typename T_Tag>
  struct TagTapeReverse;

  /// Traits for everything that can be a CoDiPack tape usually the template argument of codi::ActiveType.
  /// Possible types are codi::JacobianLinearTape, codi::ForwardEvaluation, codi::PrimalValueReuseTape, etc..
  namespace TapeTraits {

    /*******************************************************************************/
    /// @name Detection of specific real value types
    /// @{

    /// If the tape inherits from ForwardEvaluation.
    template<typename Tape, typename = void>
    struct IsForwardTape : std::false_type {};

#ifndef DOXYGEN_DISABLE
    template<typename Tape>
    struct IsForwardTape<
        Tape, typename enable_if_base_of<ForwardEvaluation<typename Tape::Real, typename Tape::Gradient>, Tape>::type>
        : std::true_type {};

    template<typename Tape>
    struct IsForwardTape<
        Tape, typename enable_if_base_of<TagTapeForward<typename Tape::Real, typename Tape::Tag>, Tape>::type>
        : std::true_type {};
#endif

#if CODI_IS_CPP14
    /// Value entry of IsForwardTape
    template<typename Tape>
    bool constexpr isForwardTape = IsForwardTape<Tape>::value;
#endif

    /// Enable if wrapper for IsForwardTape
    template<typename Tape>
    using EnableIfForwardTape = typename std::enable_if<IsForwardTape<Tape>::value>::type;

    /// If the tape inherits from PrimalValueBaseTape.
    template<typename Tape, typename = void>
    struct IsPrimalValueTape : std::false_type {};

#ifndef DOXYGEN_DISABLE
    template<typename Tape>
    struct IsPrimalValueTape<
        Tape, typename enable_if_base_of<PrimalValueBaseTape<typename Tape::TapeTypes, Tape>, Tape>::type>
        : std::true_type {};
#endif

#if CODI_IS_CPP14
    /// Value entry of IsPrimalValueTape
    template<typename Tape>
    bool constexpr isPrimalValueTape = IsPrimalValueTape<Tape>::value;
#endif

    /// Enable if wrapper for IsPrimalValueTape
    template<typename Tape>
    using EnableIfPrimalValueTape = typename std::enable_if<IsPrimalValueTape<Tape>::value>::type;

    /// If the tape inherits from JacobianBaseTape.
    template<typename Tape, typename = void>
    struct IsJacobianTape : std::false_type {};

#ifndef DOXYGEN_DISABLE
    template<typename Tape>
    struct IsJacobianTape<Tape,
                          typename enable_if_base_of<JacobianBaseTape<typename Tape::TapeTypes, Tape>, Tape>::type>
        : std::true_type {};
#endif

#if CODI_IS_CPP14
    /// Value entry of IsJacobianTape
    template<typename Tape>
    bool constexpr isJacobianTape = IsJacobianTape<Tape>::value;
#endif

    /// Enable if wrapper for IsJacobianTape
    template<typename Tape>
    using EnableIfJacobianTape = typename std::enable_if<IsJacobianTape<Tape>::value>::type;

    template<typename Tape, typename = void>
    struct IsReverseTape : std::false_type {};

#ifndef DOXYGEN_DISABLE
    /// If the tape inherits either from JacobianBaseTape or PrimalValueBaseTape.
    template<typename Tape>
    struct IsReverseTape<Tape,
                         typename std::enable_if<IsJacobianTape<Tape>::value || IsPrimalValueTape<Tape>::value>::type>
        : std::true_type {};

    template<typename Tape>
    struct IsReverseTape<
        Tape, typename enable_if_base_of<TagTapeReverse<typename Tape::Real, typename Tape::Tag>, Tape>::type>
        : std::true_type {};
#endif

#if CODI_IS_CPP14
    /// Value entry of IsReverseTape
    template<typename Tape>
    bool constexpr isReverseTape = IsReverseTape<Tape>::value;
#endif

    /// Enable if wrapper for IsReverseTape
    template<typename Tape>
    using EnableIfReverseTape = typename std::enable_if<IsReverseTape<Tape>::value>::type;

    template<typename Tape, typename = void>
    struct SupportsEditing : std::false_type {};

#ifndef DOXYGEN_DISABLE
    template<typename Tape>
    struct SupportsEditing<Tape, typename enable_if_base_of<EditingTapeInterface<typename Tape::Position>, Tape>::type>
        : std::true_type {};
#endif

#if CODI_IS_CPP14
    /// Value entry of SupportsEditing
    template<typename Tape>
    bool constexpr supportsEditing = SupportsEditing<Tape>::value;
#endif

    /// Enable if wrapper for SupportsEditing
    template<typename Tape>
    using EnableIfSupportsEditing = typename std::enable_if<SupportsEditing<Tape>::value>::type;

    /// Enable if wrapper for SupportsEditing
    template<typename Tape>
    using EnableIfNoEditing = typename std::enable_if<!SupportsEditing<Tape>::value>::type;

    /// @}
  }
}
