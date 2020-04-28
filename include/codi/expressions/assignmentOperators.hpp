#pragma once

#include "../config.h"
#include "../aux/macros.h"
#include "../tapes/interfaces/internalExpressionTapeInterface.hpp"
#include "lhsExpressionInterface.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  template<typename _Tape, typename _Impl>
  struct AssignmentOperators {
    public:

      using Tape = DECLARE_DEFAULT(_Tape, TEMPLATE(InternalExpressionTapeInterface<int>));
      using Impl = DECLARE_DEFAULT(_Impl, TEMPLATE(LhsExpressionInterface<double, int, Tape, _Impl>));

      using Real = DECLARE_DEFAULT(typename Tape::Real, double);
      using PassiveReal = PassiveRealType<Real>;

      CODI_INLINE Impl& cast() {
        return static_cast<Impl&>(*this);
      }

      template<typename Rhs>
      CODI_INLINE Impl& operator+=(ExpressionInterface<Real, Rhs> const& rhs) {
        return cast() = (cast() + rhs);
      }

      template<typename Rhs>
      CODI_INLINE Impl& operator-=(ExpressionInterface<Real, Rhs> const& rhs) {
        return cast() = (cast() - rhs);
      }

      template<typename Rhs>
      CODI_INLINE Impl& operator*=(ExpressionInterface<Real, Rhs> const& rhs) {
        return cast() = (cast() * rhs);
      }

      template<typename Rhs>
      CODI_INLINE Impl& operator/=(ExpressionInterface<Real, Rhs> const& rhs) {
        return cast() = (cast() / rhs);
      }

      CODI_INLINE Impl& operator+=(PassiveReal const& rhs) {
        if(Tape::AllowJacobianOptimization) {
          cast().value() += rhs;
        } else {
          cast() = (cast() + rhs);
        }
        return cast();
      }

      CODI_INLINE Impl& operator-=(PassiveReal const& rhs) {
        if(Tape::AllowJacobianOptimization) {
          cast().value() -= rhs;
        } else {
          cast() = (cast() - rhs);
        }
        return cast();
      }

      CODI_INLINE Impl& operator*=(PassiveReal const& rhs) {
        return cast() = (cast() * rhs);
      }

      CODI_INLINE Impl& operator/=(PassiveReal const& rhs) {
        return cast() = (cast() / rhs);
      }
  };
}
