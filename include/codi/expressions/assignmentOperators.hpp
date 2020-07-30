#pragma once

#include "../config.h"
#include "../aux/macros.hpp"
#include "../tapes/interfaces/internalExpressionTapeInterface.hpp"
#include "lhsExpressionInterface.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  template<typename _Tape, typename _Impl>
  struct AssignmentOperators {
    public:

      using Tape = CODI_DECLARE_DEFAULT(_Tape, CODI_TEMPLATE(InternalExpressionTapeInterface<int>));
      using Impl = CODI_DECLARE_DEFAULT(_Impl, CODI_TEMPLATE(LhsExpressionInterface<double, int, Tape, _Impl>));

      using Real = CODI_DECLARE_DEFAULT(typename Tape::Real, double);
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
