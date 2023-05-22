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

#include "../config.h"
#include "../misc/macros.hpp"
#include "../tapes/interfaces/internalStatementRecordingTapeInterface.hpp"
#include "lhsExpressionInterface.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Implementation of assignment operators for LhsExpressionInterface implementations.
   *
   * Implements: +=, -=, *=, /= for Expressions and passive values.
   *
   * @tparam T_Tape  The tape of the lvalue implementation.
   * @tparam T_Impl  The lvalue LhsExpressionInterface implementation.
   */
  template<typename T_Tape, typename T_Impl>
  struct AssignmentOperators {
    public:

      using Tape = CODI_DD(T_Tape, CODI_T(InternalStatementRecordingTapeInterface<int>));  ///< See AssignmentOperators.
      using Impl = CODI_DD(T_Impl,
                           CODI_T(LhsExpressionInterface<double, double, Tape, T_Impl>));  ///< See AssignmentOperators.

      using Real = CODI_DD(typename Tape::Real, double);  ///< See InternalStatementRecordingTapeInterface.
      using PassiveReal = RealTraits::PassiveReal<Real>;  ///< Basic computation type.

      /// Cast to the implementation.
      CODI_INLINE Impl& cast() {
        return static_cast<Impl&>(*this);
      }

      /// Operator += for expressions.
      template<typename Rhs>
      CODI_INLINE Impl& operator+=(ExpressionInterface<Real, Rhs> const& rhs) {
        return cast() = (cast() + rhs);
      }

      /// Operator -= for expressions.
      template<typename Rhs>
      CODI_INLINE Impl& operator-=(ExpressionInterface<Real, Rhs> const& rhs) {
        return cast() = (cast() - rhs);
      }

      /// Operator *= for expressions.
      template<typename Rhs>
      CODI_INLINE Impl& operator*=(ExpressionInterface<Real, Rhs> const& rhs) {
        return cast() = (cast() * rhs);
      }

      /// Operator /= for expressions.
      template<typename Rhs>
      CODI_INLINE Impl& operator/=(ExpressionInterface<Real, Rhs> const& rhs) {
        return cast() = (cast() / rhs);
      }

      /// Operator += for passive values.
      CODI_INLINE Impl& operator+=(PassiveReal const& rhs) {
        if (Tape::AllowJacobianOptimization) {
          cast().value() += rhs;
        } else {
          cast() = (cast() + rhs);
        }
        return cast();
      }

      /// Operator -= for passive values.
      CODI_INLINE Impl& operator-=(PassiveReal const& rhs) {
        if (Tape::AllowJacobianOptimization) {
          cast().value() -= rhs;
        } else {
          cast() = (cast() - rhs);
        }
        return cast();
      }

      /// Operator *= for passive values.
      CODI_INLINE Impl& operator*=(PassiveReal const& rhs) {
        return cast() = (cast() * rhs);
      }

      /// Operator /= for passive values.
      CODI_INLINE Impl& operator/=(PassiveReal const& rhs) {
        return cast() = (cast() / rhs);
      }
  };
}
