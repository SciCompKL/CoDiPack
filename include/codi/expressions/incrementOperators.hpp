#pragma once

#include "../config.h"
#include "../aux/macros.h"
#include "../tapes/interfaces/internalStatementRecordingInterface.hpp"
#include "lhsExpressionInterface.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  template<typename _Tape, typename _Impl>
  struct IncrementOperators {
    public:

      using Tape = DECLARE_DEFAULT(_Tape, TEMPLATE(InternalStatementRecordingInterface<int>));
      using Impl = DECLARE_DEFAULT(_Impl, TEMPLATE(LhsExpressionInterface<double, int, Tape, _Impl>));

      using Real = typename Tape::Real;
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
