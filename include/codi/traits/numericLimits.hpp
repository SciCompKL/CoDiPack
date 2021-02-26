#pragma once

#include <limits>

#include "../aux/macros.hpp"
#include "../config.h"
#include "realTraits.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  template<typename _Tape>
  struct ActiveType;
}

namespace std {

  /// Specialization of std::numeric_limits for CoDiPack active types.
  template<typename Tape>
  struct numeric_limits<codi::ActiveType<Tape>> {
    private:

      using Type = codi::ActiveType<Tape>;
      using Passive = codi::RealTraits::PassiveReal<Type>;

    public:

      static bool constexpr is_specialized = true;  ///< See numeric_limits
      static Type constexpr min() {
        return Type(numeric_limits<Passive>::min());
      }  ///< See numeric_limits
      static Type constexpr max() {
        return Type(numeric_limits<Passive>::max());
      }  ///< See numeric_limits
      static Type constexpr lowest() {
        return Type(numeric_limits<Passive>::lowest());
      }  ///< See numeric_limits
      static int constexpr digits = numeric_limits<Passive>::digits;  ///< See numeric_limits
      static int constexpr digits10 = numeric_limits<Passive>::digits10;  ///< See numeric_limits
      static bool constexpr is_signed = numeric_limits<Passive>::is_signed;  ///< See numeric_limits
      static bool constexpr is_integer = numeric_limits<Passive>::is_integer;  ///< See numeric_limits
      static bool constexpr is_exact = numeric_limits<Passive>::is_exact;  ///< See numeric_limits
      static int constexpr radix = numeric_limits<Passive>::is_exact;  ///< See numeric_limits
      static Type constexpr epsilon() {
        return Type(numeric_limits<Passive>::epsilon());
      }  ///< See numeric_limits
      static Type constexpr round_error() {
        return Type(numeric_limits<Passive>::round_error());
      }  ///< See numeric_limits

      static int constexpr min_exponent = numeric_limits<Passive>::min_exponent;  ///< See numeric_limits
      static int constexpr min_exponent10 = numeric_limits<Passive>::max_exponent10;  ///< See numeric_limits
      static int constexpr max_exponent = numeric_limits<Passive>::max_exponent;  ///< See numeric_limits
      static int constexpr max_exponent10 = numeric_limits<Passive>::max_exponent10;  ///< See numeric_limits

      static bool constexpr has_infinity = numeric_limits<Passive>::has_infinity;  ///< See numeric_limits
      static bool constexpr has_quiet_NaN = numeric_limits<Passive>::has_quiet_NaN;  ///< See numeric_limits
      static bool constexpr has_signaling_NaN = numeric_limits<Passive>::has_signaling_NaN;  ///< See numeric_limits
      static float_denorm_style constexpr has_denorm = numeric_limits<Passive>::has_denorm;  ///< See numeric_limits
      static bool constexpr has_denorm_loss = numeric_limits<Passive>::has_denorm_loss;  ///< See numeric_limits
      static Type constexpr infinity() {
        return Type(numeric_limits<Passive>::infinity());
      }  ///< See numeric_limits
      static Type constexpr quiet_NaN() {
        return Type(numeric_limits<Passive>::quiet_NaN());
      }  ///< See numeric_limits
      static Type constexpr signaling_NaN() {
        return Type(numeric_limits<Passive>::signaling_NaN());
      }  ///< See numeric_limits
      static Type constexpr denorm_min() {
        return Type(numeric_limits<Passive>::denorm_min());
      }  ///< See numeric_limits

      static bool constexpr is_iec559 = numeric_limits<Passive>::is_iec559;  ///< See numeric_limits
      static bool constexpr is_bounded = numeric_limits<Passive>::is_bounded;  ///< See numeric_limits
      static bool constexpr is_modulo = numeric_limits<Passive>::is_modulo;  ///< See numeric_limits

      static bool constexpr traps = numeric_limits<Passive>::traps;  ///< See numeric_limits
      static bool constexpr tinyness_before = numeric_limits<Passive>::tinyness_before;  ///< See numeric_limits
      static float_round_style constexpr round_style = numeric_limits<Passive>::round_style;  ///< See numeric_limits
  };
}
