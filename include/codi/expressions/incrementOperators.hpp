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
   * Implements: prefix ++, postfix ++, prefix --, postfix --
   *
   * @tparam T_Tape  The tape of the lvalue implementation.
   * @tparam T_Impl  The lvalue LhsExpressionInterface implementation.
   */
  template<typename T_Tape, typename T_Impl>
  struct IncrementOperators {
    public:

      using Tape = CODI_DD(T_Tape, CODI_T(InternalStatementRecordingTapeInterface<int>));  ///< See IncrementOperators.
      using Impl = CODI_DD(T_Impl,
                           CODI_T(LhsExpressionInterface<double, int, Tape, T_Impl>));  ///< See IncrementOperators.

      using Real = CODI_DD(typename Tape::Real, double);  ///< See InternalStatementRecordingTapeInterface.
      using PassiveReal = RealTraits::PassiveReal<Real>;  ///< Basic computation type.

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
