/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2026 Chair for Scientific Computing (SciComp), RPTU University Kaiserslautern-Landau
 * Homepage: http://scicomp.rptu.de
 * Contact:  Prof. Nicolas R. Gauger (codi@scicomp.uni-kl.de)
 *
 * Lead developers: Max Sagebaum, Johannes Blühdorn (SciComp, RPTU University Kaiserslautern-Landau)
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
 *  - SciComp, RPTU University Kaiserslautern-Landau:
 *    - Max Sagebaum
 *    - Johannes Blühdorn
 *    - Former members:
 *      - Tim Albring
 */
#pragma once

#include "../config.h"
#include "../misc/macros.hpp"
#include "lhsExpressionInterface.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /// Default implementations for the passive += and -= operators.
  template<typename T_Real, bool T_JacobianOptimization, typename T_Impl>
  struct AssignmentOperatorsPassiveJacobianOptimization {
    public:

      using Real = CODI_DD(T_Real, double);  ///< See AssignmentOperators.
      static bool constexpr JacobianOptimization =
          CODI_DD(T_JacobianOptimization, false);  ///< See AssignmentOperators.
      using Impl = CODI_DD(T_Impl,
                           CODI_T(LhsExpressionInterface<double, int, Tape, T_Impl>));  ///< See AssignmentOperators.

      using PassiveReal = RealTraits::PassiveReal<Real>;  ///< Basic computation type.

      /// Operator += for passive values.
      CODI_INLINE Impl& operator+=(PassiveReal const& rhs) {
        cast() = (cast() + rhs);
        return cast();
      }

      /// Operator -= for passive values.
      CODI_INLINE Impl& operator-=(PassiveReal const& rhs) {
        cast() = (cast() - rhs);
        return cast();
      }

    private:

      /// Cast to the implementation.
      CODI_INLINE Impl& cast() {
        return static_cast<Impl&>(*this);
      }
  };

  /// Optimized Jacobian implementations for the passive += and -= operators.
  template<typename T_Real, typename T_Impl>
  struct AssignmentOperatorsPassiveJacobianOptimization<T_Real, true, T_Impl> {
    public:

      using Real = CODI_DD(T_Real, double);               ///< See AssignmentOperators.
      static bool constexpr JacobianOptimization = true;  ///< See AssignmentOperators.
      using Impl = CODI_DD(T_Impl,
                           CODI_T(LhsExpressionInterface<double, int, Tape, T_Impl>));  ///< See AssignmentOperators.

      using PassiveReal = RealTraits::PassiveReal<Real>;  ///< Basic computation type.

      /// Operator += for passive values.
      CODI_INLINE Impl& operator+=(PassiveReal const& rhs) {
        cast().value() += rhs;
        return cast();
      }

      /// Operator -= for passive values.
      CODI_INLINE Impl& operator-=(PassiveReal const& rhs) {
        cast().value() -= rhs;
        return cast();
      }

    private:

      /// Cast to the implementation.
      CODI_INLINE Impl& cast() {
        return static_cast<Impl&>(*this);
      }
  };

  /**
   * @brief Provides assignment operators for LhsExpressionInterface implementations.
   *
   * Implements: +=, -=, *=, /= for Expressions and passive values.
   *
   * @tparam T_Real  The real type for the right hand side expressions.
   * @tparam T_JacobianOptimization  If Jacobian optimization is allows for += and -= operators.
   * @tparam T_Impl  The lvalue LhsExpressionInterface implementation.
   */
  template<typename T_Real, bool T_JacobianOptimization, typename T_Impl>
  struct AssignmentOperators
      : public AssignmentOperatorsPassiveJacobianOptimization<T_Real, T_JacobianOptimization, T_Impl> {
    public:

      using Real = CODI_DD(T_Real, double);                                                ///< See AssignmentOperators.
      static bool constexpr JacobianOptimization = CODI_DD(T_JacobianOptimization, true);  ///< See AssignmentOperators.
      using Impl = CODI_DD(T_Impl,
                           CODI_T(LhsExpressionInterface<double, double, Tape, T_Impl>));  ///< See AssignmentOperators.

      using Base =
          AssignmentOperatorsPassiveJacobianOptimization<T_Real, T_JacobianOptimization, T_Impl>;  ///< Abbreviation for
                                                                                                   ///< base class.
      using PassiveReal = RealTraits::PassiveReal<Real>;  ///< Basic computation type.

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

      using Base::operator+=;
      using Base::operator-=;

      /// Operator *= for passive values.
      CODI_INLINE Impl& operator*=(PassiveReal const& rhs) {
        return cast() = (cast() * rhs);
      }

      /// Operator /= for passive values.
      CODI_INLINE Impl& operator/=(PassiveReal const& rhs) {
        return cast() = (cast() / rhs);
      }

    private:

      /// Cast to the implementation.
      CODI_INLINE Impl& cast() {
        return static_cast<Impl&>(*this);
      }
  };
}
