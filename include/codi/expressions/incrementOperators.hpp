#pragma once

#include "../config.h"
#include "../aux/macros.hpp"
#include "../tapes/interfaces/internalStatementRecordingInterface.hpp"
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

      using Tape = CODI_DECLARE_DEFAULT(_Tape, CODI_TEMPLATE(InternalStatementRecordingInterface<int>)); ///< See IncrementOperators
      using Impl = CODI_DECLARE_DEFAULT(_Impl, CODI_TEMPLATE(LhsExpressionInterface<double, int, Tape, _Impl>)); ///< See IncrementOperators

      using Real = CODI_DECLARE_DEFAULT(typename Tape::Real, double); ///< See InternalStatementRecordingInterface
      using PassiveReal = PassiveRealType<Real>; ///< Basic computation type

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
