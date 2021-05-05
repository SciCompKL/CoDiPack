#pragma once

#include "../aux/macros.hpp"
#include "../config.h"
#include "../tapes/interfaces/internalStatementRecordingTapeInterface.hpp"
#include "lhsExpressionInterface.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Implementation of increment operators for LhsExpressionInterface implementations.
   *
   * Implements: prefix ++, postfix, ++ prefix --, postfix --
   *
   * @tparam _Tape  The tape of the lvalue implementation.
   * @tparam _Impl  The lvalue LhsExpressionInterface implementation.
   */
  template<typename _Tape, typename _Impl>
  struct IncrementOperators {
    public:

      using Tape = CODI_DD(_Tape, CODI_T(InternalStatementRecordingTapeInterface<int>));  ///< See IncrementOperators
      using Impl = CODI_DD(_Impl,
                           CODI_T(LhsExpressionInterface<double, int, Tape, _Impl>));  ///< See IncrementOperators

      using Real = CODI_DD(typename Tape::Real, double);  ///< See InternalStatementRecordingTapeInterface
      using PassiveReal = RealTraits::PassiveReal<Real>;  ///< Basic computation type

      /// Cast to the implementation.
      CODI_INLINE Impl& cast() {
        return static_cast<Impl&>(*this);
      }

      /// Prefix operator ++
      CODI_INLINE Impl& operator++() {
        return cast() = cast() + PassiveReal(1.0);
      }

      /// Postfix operator ++
      CODI_INLINE Impl operator++(int u) {
        CODI_UNUSED(u);

        Impl r(cast());
        cast() = cast() + PassiveReal(1.0);
        return r;
      }

      /// Prefix operator --
      CODI_INLINE Impl& operator--() {
        return cast() = cast() - PassiveReal(1.0);
      }

      /// Postfix operator --
      CODI_INLINE Impl operator--(int u) {
        CODI_UNUSED(u);

        Impl r(cast());
        cast() = cast() - PassiveReal(1.0);
        return r;
      }
  };
}
