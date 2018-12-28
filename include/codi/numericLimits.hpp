/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2019 Chair for Scientific Computing (SciComp), TU Kaiserslautern
 * Homepage: http://www.scicomp.uni-kl.de
 * Contact:  Prof. Nicolas R. Gauger (codi@scicomp.uni-kl.de)
 *
 * Lead developers: Max Sagebaum, Tim Albring (SciComp, TU Kaiserslautern)
 *
 * This file is part of CoDiPack (http://www.scicomp.uni-kl.de/software/codi).
 *
 * CoDiPack is free software: you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * CoDiPack is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty
 * of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 * See the GNU General Public License for more details.
 * You should have received a copy of the GNU
 * General Public License along with CoDiPack.
 * If not, see <http://www.gnu.org/licenses/>.
 *
 * Authors: Max Sagebaum, Tim Albring, (SciComp, TU Kaiserslautern)
 */

#pragma once

#include <limits>

#include "activeReal.hpp"
#include "typeTraits.hpp"

namespace std {

  /**
   * @brief The numeric limits for the ActiveReal's are based on the numeric limits for the
   *        passive type (aka the start of the ActiveReal chain).
   *
   * @tparam Tape  The tape that defines the properties of the ActiveReal.
   */
  template <typename Tape>
  class numeric_limits<codi::ActiveReal<Tape> > {
    private:
      /** @brief The active type for which the numeric limits are specialized */
      typedef codi::ActiveReal<Tape> Real;

      /** @brief The passive type from which the values are taken */
      typedef typename codi::TypeTraits<codi::ActiveReal<Tape> >::PassiveReal Passive;
    public:

      /** @brief Needs to be false per definition */
      static constexpr bool is_specialized = false;
      /** @brief Use the value from the passive type */
      static constexpr Real min() { return Real(numeric_limits<Passive>::min()); }
      /** @brief Use the value from the passive type */
      static constexpr Real max() { return Real(numeric_limits<Passive>::max()); }
      /** @brief Use the value from the passive type */
      static constexpr Real lowest() { return Real(numeric_limits<Passive>::lowest()); }
      /** @brief Use the value from the passive type */
      static constexpr int  digits = numeric_limits<Passive>::digits;
      /** @brief Use the value from the passive type */
      static constexpr int  digits10 = numeric_limits<Passive>::digits10;
      /** @brief Use the value from the passive type */
      static constexpr bool is_signed = numeric_limits<Passive>::is_signed;
      /** @brief Use the value from the passive type */
      static constexpr bool is_integer = numeric_limits<Passive>::is_integer;
      /** @brief Use the value from the passive type */
      static constexpr bool is_exact = numeric_limits<Passive>::is_exact;
      /** @brief Use the value from the passive type */
      static constexpr int radix = numeric_limits<Passive>::is_exact;
      /** @brief Use the value from the passive type */
      static constexpr Real epsilon() { return Real(numeric_limits<Passive>::epsilon()); }
      /** @brief Use the value from the passive type */
      static constexpr Real round_error() { return Real(numeric_limits<Passive>::round_error()); }

      /** @brief Use the value from the passive type */
      static constexpr int  min_exponent = numeric_limits<Passive>::min_exponent;
      /** @brief Use the value from the passive type */
      static constexpr int  min_exponent10 = numeric_limits<Passive>::max_exponent10;
      /** @brief Use the value from the passive type */
      static constexpr int  max_exponent = numeric_limits<Passive>::max_exponent;
      /** @brief Use the value from the passive type */
      static constexpr int  max_exponent10 = numeric_limits<Passive>::max_exponent10;

      /** @brief Use the value from the passive type */
      static constexpr bool has_infinity = numeric_limits<Passive>::has_infinity;
      /** @brief Use the value from the passive type */
      static constexpr bool has_quiet_NaN = numeric_limits<Passive>::has_quiet_NaN;
      /** @brief Use the value from the passive type */
      static constexpr bool has_signaling_NaN = numeric_limits<Passive>::has_signaling_NaN;
      /** @brief Use the value from the passive type */
      static constexpr float_denorm_style has_denorm = numeric_limits<Passive>::has_denorm;
      /** @brief Use the value from the passive type */
      static constexpr bool has_denorm_loss = numeric_limits<Passive>::has_denorm_loss;
      /** @brief Use the value from the passive type */
      static constexpr Real infinity() { return Real(numeric_limits<Passive>::infinity()); }
      /** @brief Use the value from the passive type */
      static constexpr Real quiet_NaN() { return Real(numeric_limits<Passive>::quiet_NaN()); }
      /** @brief Use the value from the passive type */
      static constexpr Real signaling_NaN() { return Real(numeric_limits<Passive>::signaling_NaN()); }
      /** @brief Use the value from the passive type */
      static constexpr Real denorm_min() { return Real(numeric_limits<Passive>::denorm_min()); }

      /** @brief Use the value from the passive type */
      static constexpr bool is_iec559 = numeric_limits<Passive>::is_iec559;
      /** @brief Use the value from the passive type */
      static constexpr bool is_bounded = numeric_limits<Passive>::is_bounded;
      /** @brief Use the value from the passive type */
      static constexpr bool is_modulo = numeric_limits<Passive>::is_modulo;

      /** @brief Use the value from the passive type */
      static constexpr bool traps = numeric_limits<Passive>::traps;
      /** @brief Use the value from the passive type */
      static constexpr bool tinyness_before = numeric_limits<Passive>::tinyness_before;
      /** @brief Use the value from the passive type */
      static constexpr float_round_style round_style = numeric_limits<Passive>::round_style;
  };
}
