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

      static constexpr bool is_specialized = true;
      static constexpr Type min() { return Type(numeric_limits<Passive>::min()); }
      static constexpr Type max() { return Type(numeric_limits<Passive>::max()); }
      static constexpr Type lowest() { return Type(numeric_limits<Passive>::lowest()); }
      static constexpr int  digits = numeric_limits<Passive>::digits;
      static constexpr int  digits10 = numeric_limits<Passive>::digits10;
      static constexpr bool is_signed = numeric_limits<Passive>::is_signed;
      static constexpr bool is_integer = numeric_limits<Passive>::is_integer;
      static constexpr bool is_exact = numeric_limits<Passive>::is_exact;
      static constexpr int radix = numeric_limits<Passive>::is_exact;
      static constexpr Type epsilon() { return Type(numeric_limits<Passive>::epsilon()); }
      static constexpr Type round_error() { return Type(numeric_limits<Passive>::round_error()); }

      static constexpr int  min_exponent = numeric_limits<Passive>::min_exponent;
      static constexpr int  min_exponent10 = numeric_limits<Passive>::max_exponent10;
      static constexpr int  max_exponent = numeric_limits<Passive>::max_exponent;
      static constexpr int  max_exponent10 = numeric_limits<Passive>::max_exponent10;

      static constexpr bool has_infinity = numeric_limits<Passive>::has_infinity;
      static constexpr bool has_quiet_NaN = numeric_limits<Passive>::has_quiet_NaN;
      static constexpr bool has_signaling_NaN = numeric_limits<Passive>::has_signaling_NaN;
      static constexpr float_denorm_style has_denorm = numeric_limits<Passive>::has_denorm;
      static constexpr bool has_denorm_loss = numeric_limits<Passive>::has_denorm_loss;
      static constexpr Type infinity() { return Type(numeric_limits<Passive>::infinity()); }
      static constexpr Type quiet_NaN() { return Type(numeric_limits<Passive>::quiet_NaN()); }
      static constexpr Type signaling_NaN() { return Type(numeric_limits<Passive>::signaling_NaN()); }
      static constexpr Type denorm_min() { return Type(numeric_limits<Passive>::denorm_min()); }

      static constexpr bool is_iec559 = numeric_limits<Passive>::is_iec559;
      static constexpr bool is_bounded = numeric_limits<Passive>::is_bounded;
      static constexpr bool is_modulo = numeric_limits<Passive>::is_modulo;

      static constexpr bool traps = numeric_limits<Passive>::traps;
      static constexpr bool tinyness_before = numeric_limits<Passive>::tinyness_before;
      static constexpr float_round_style round_style = numeric_limits<Passive>::round_style;
  };
}
