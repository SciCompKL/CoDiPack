#pragma once

#include "../config.h"
#include "../aux/macros.hpp"
#include "../tapes/interfaces/internalStatementRecordingInterface.hpp"
#include "lhsExpressionInterface.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Implementation of assignment operators for LhsExpressionInterface implementations.
   *
   * Implements: +=, -=, *=, /= for Expressions and passive values.
   *
   * @tparam _Tape  The tape of the lvalue implementation.
   * @tparam _Impl  The lvalue LhsExpressionInterface implementation.
   */
  template<typename _Tape, typename _Impl>
  struct AssignmentOperators {
    public:

      using Tape = CODI_DECLARE_DEFAULT(_Tape, CODI_TEMPLATE(InternalStatementRecordingInterface<int>)); ///< See AssignmentOperators
      using Impl = CODI_DECLARE_DEFAULT(_Impl, CODI_TEMPLATE(LhsExpressionInterface<double, int, Tape, _Impl>)); ///< See AssignmentOperators

      using Real = CODI_DECLARE_DEFAULT(typename Tape::Real, double); ///< See InternalStatementRecordingInterface
      using PassiveReal = RealTraits::PassiveReal<Real>; ///< Basic computation type

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
        if(Tape::AllowJacobianOptimization) {
          cast().value() += rhs;
        } else {
          cast() = (cast() + rhs);
        }
        return cast();
      }

      /// Operator -= for passive values.
      CODI_INLINE Impl& operator-=(PassiveReal const& rhs) {
        if(Tape::AllowJacobianOptimization) {
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
