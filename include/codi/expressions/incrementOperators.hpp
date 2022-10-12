/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2022 Chair for Scientific Computing (SciComp), TU Kaiserslautern
 * Homepage: http://www.scicomp.uni-kl.de
 * Contact:  Prof. Nicolas R. Gauger (codi@scicomp.uni-kl.de)
 *
 * Lead developers: Max Sagebaum, Johannes Blühdorn (SciComp, TU Kaiserslautern)
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
 *  - SciComp, TU Kaiserslautern:
 *    - Max Sagebaum
 *    - Johannes Blühdorn
 *    - Former members:
 *      - Tim Albring
 */
#pragma once

#include "../misc/macros.hpp"
#include "../config.h"
#include "../tapes/interfaces/internalStatementRecordingTapeInterface.hpp"
#include "lhsExpressionInterface.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Implementation of increment operators for LhsExpressionInterface implementations.
   *
   * Implements: prefix ++, postfix ++, prefix --, postfix --
   *
   * @tparam T_Tape  The tape of the lvalue implementation.
   * @tparam T_Impl  The lvalue LhsExpressionInterface implementation.
   */
  template<typename T_Tape, typename T_Impl>
  struct IncrementOperators {
    public:

      using Tape = CODI_DD(T_Tape, CODI_T(InternalStatementRecordingTapeInterface<int>));  ///< See IncrementOperators.
      using Impl = CODI_DD(T_Impl, CODI_T(LhsExpressionInterface<double, double,
                                                                 Tape, T_Impl>));  ///< See IncrementOperators.

      using Real = CODI_DD(typename Tape::Real, double);  ///< See InternalStatementRecordingTapeInterface.
      using PassiveReal = RealTraits::PassiveReal<Real>;  ///< Basic computation type.

      /// Cast to the implementation.
      CODI_INLINE Impl& cast() {
        return static_cast<Impl&>(*this);
      }

      /// Prefix operator++
      CODI_INLINE Impl& operator++() {
        return cast() = cast() + PassiveReal(1.0);
      }

      /// Postfix operator++
      CODI_INLINE Impl operator++(int u) {
        CODI_UNUSED(u);

        Impl r(cast());
        cast() = cast() + PassiveReal(1.0);
        return r;
      }

      /// Prefix operator--
      CODI_INLINE Impl& operator--() {
        return cast() = cast() - PassiveReal(1.0);
      }

      /// Postfix operator--
      CODI_INLINE Impl operator--(int u) {
        CODI_UNUSED(u);

        Impl r(cast());
        cast() = cast() - PassiveReal(1.0);
        return r;
      }
  };
}
