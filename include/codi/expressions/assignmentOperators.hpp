#pragma once

#include "../aux/macros.hpp"
#include "../config.h"
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
                           CODI_T(LhsExpressionInterface<double, int, Tape, T_Impl>));  ///< See AssignmentOperators.

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
