#pragma once

#include <limits>

#include "../aux/macros.h"
#include "../config.h"
#include "realTraits.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  template<typename _Tape>
  struct ActiveType;
}

namespace std {

  template <typename Tape>
  struct numeric_limits<codi::ActiveType<Tape>> {
    private:

      using Type = codi::ActiveType<Tape>;
      using Passive = codi::PassiveRealType<Type>;

    public:

      static bool constexpr is_specialized = true;
      static Type constexpr min() { return Type(numeric_limits<Passive>::min()); }
      static Type constexpr max() { return Type(numeric_limits<Passive>::max()); }
      static Type constexpr lowest() { return Type(numeric_limits<Passive>::lowest()); }
      static int  constexpr digits = numeric_limits<Passive>::digits;
      static int  constexpr digits10 = numeric_limits<Passive>::digits10;
      static bool constexpr is_signed = numeric_limits<Passive>::is_signed;
      static bool constexpr is_integer = numeric_limits<Passive>::is_integer;
      static bool constexpr is_exact = numeric_limits<Passive>::is_exact;
      static int  constexpr radix = numeric_limits<Passive>::is_exact;
      static Type constexpr epsilon() { return Type(numeric_limits<Passive>::epsilon()); }
      static Type constexpr round_error() { return Type(numeric_limits<Passive>::round_error()); }

      static int  constexpr min_exponent = numeric_limits<Passive>::min_exponent;
      static int  constexpr min_exponent10 = numeric_limits<Passive>::max_exponent10;
      static int  constexpr max_exponent = numeric_limits<Passive>::max_exponent;
      static int  constexpr max_exponent10 = numeric_limits<Passive>::max_exponent10;

      static bool constexpr has_infinity = numeric_limits<Passive>::has_infinity;
      static bool constexpr has_quiet_NaN = numeric_limits<Passive>::has_quiet_NaN;
      static bool constexpr has_signaling_NaN = numeric_limits<Passive>::has_signaling_NaN;
      static float_denorm_style constexpr has_denorm = numeric_limits<Passive>::has_denorm;
      static bool constexpr has_denorm_loss = numeric_limits<Passive>::has_denorm_loss;
      static Type constexpr infinity() { return Type(numeric_limits<Passive>::infinity()); }
      static Type constexpr quiet_NaN() { return Type(numeric_limits<Passive>::quiet_NaN()); }
      static Type constexpr signaling_NaN() { return Type(numeric_limits<Passive>::signaling_NaN()); }
      static Type constexpr denorm_min() { return Type(numeric_limits<Passive>::denorm_min()); }

      static bool constexpr is_iec559 = numeric_limits<Passive>::is_iec559;
      static bool constexpr is_bounded = numeric_limits<Passive>::is_bounded;
      static bool constexpr is_modulo = numeric_limits<Passive>::is_modulo;

      static bool constexpr traps = numeric_limits<Passive>::traps;
      static bool constexpr tinyness_before = numeric_limits<Passive>::tinyness_before;
      static float_round_style constexpr round_style = numeric_limits<Passive>::round_style;
  };
}
