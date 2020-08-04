#pragma once

#include "../config.h"
#include "../aux/macros.hpp"
#include "../tapes/interfaces/internalStatementRecordingInterface.hpp"
#include "lhsExpressionInterface.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  template<typename _Tape, typename _Impl>
  struct IncrementOperators {
    public:

      using Tape = CODI_DECLARE_DEFAULT(_Tape, CODI_TEMPLATE(InternalStatementRecordingInterface<int>));
      using Impl = CODI_DECLARE_DEFAULT(_Impl, CODI_TEMPLATE(LhsExpressionInterface<double, int, Tape, _Impl>));

      using Real = CODI_DECLARE_DEFAULT(typename Tape::Real, double);
      using PassiveReal = PassiveRealType<Real>;

      CODI_INLINE Impl& cast() {
        return static_cast<Impl&>(*this);
      }

      CODI_INLINE Impl& operator++() {
        return cast() = cast() + PassiveReal(1.0);
      }

      CODI_INLINE Impl operator++(int u) {
        CODI_UNUSED(u);

        Impl r(cast());
        cast() = cast() + PassiveReal(1.0);
        return r;
      }

      CODI_INLINE Impl& operator--() {
        return cast() = cast() - PassiveReal(1.0);
      }

      CODI_INLINE Impl operator--(int u) {
        CODI_UNUSED(u);

        Impl r(cast());
        cast() = cast() - PassiveReal(1.0);
        return r;
      }
  };
}
