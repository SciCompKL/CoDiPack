/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015 Chair for Scientific Computing (SciComp), TU Kaiserslautern
 * Homepage: http://www.scicomp.uni-kl.de
 * Contact:  Prof. Nicolas R. Gauger (codi@scicomp.uni-kl.de)
 *
 * Lead developers: Max Sagebaum, Tim Albring (SciComp, TU Kaiserslautern)
 *
 * This file is part of CoDiPack (http://www.scicomp.uni-kl.de/software/codi).
 *
 * CoDiPack is free software: you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation, either version 2 of the
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

  template <typename Tape> class numeric_limits<codi::ActiveReal<Tape> > {

      typedef codi::ActiveReal<Tape> T;
      typedef typename codi::TypeTraits<codi::ActiveReal<Tape> >::PassiveReal Real;
  public:
    static constexpr bool is_specialized = false;
    static constexpr T min() noexcept { return T(numeric_limits<Real>::min()); }
    static constexpr T max() noexcept { return T(numeric_limits<Real>::max()); }
    static constexpr T lowest() noexcept { return T(numeric_limits<Real>::lowest()); }
    static constexpr int  digits = numeric_limits<Real>::digits;
    static constexpr int  digits10 = numeric_limits<Real>::digits10;
    static constexpr bool is_signed = numeric_limits<Real>::is_signed;
    static constexpr bool is_integer = numeric_limits<Real>::is_integer;
    static constexpr bool is_exact = numeric_limits<Real>::is_exact;
    static constexpr int radix = numeric_limits<Real>::is_exact;
    static constexpr T epsilon() noexcept { return T(numeric_limits<Real>::epsilon()); }
    static constexpr T round_error() noexcept { return T(numeric_limits<Real>::round_error()); }

    static constexpr int  min_exponent = numeric_limits<Real>::min_exponent;
    static constexpr int  min_exponent10 = numeric_limits<Real>::max_exponent10;
    static constexpr int  max_exponent = numeric_limits<Real>::max_exponent;
    static constexpr int  max_exponent10 = numeric_limits<Real>::max_exponent10;

    static constexpr bool has_infinity = numeric_limits<Real>::has_infinity;
    static constexpr bool has_quiet_NaN = numeric_limits<Real>::has_quiet_NaN;
    static constexpr bool has_signaling_NaN = numeric_limits<Real>::has_signaling_NaN;
    static constexpr float_denorm_style has_denorm = numeric_limits<Real>::has_denorm;
    static constexpr bool has_denorm_loss = numeric_limits<Real>::has_denorm_loss;
    static constexpr T infinity() noexcept { return T(numeric_limits<Real>::infinity()); }
    static constexpr T quiet_NaN() noexcept { return T(numeric_limits<Real>::quiet_NaN()); }
    static constexpr T signaling_NaN() noexcept { return T(numeric_limits<Real>::signaling_NaN()); }
    static constexpr T denorm_min() noexcept { return T(numeric_limits<Real>::denorm_min()); }

    static constexpr bool is_iec559 = numeric_limits<Real>::is_iec559;
    static constexpr bool is_bounded = numeric_limits<Real>::is_bounded;
    static constexpr bool is_modulo = numeric_limits<Real>::is_modulo;

    static constexpr bool traps = numeric_limits<Real>::traps;
    static constexpr bool tinyness_before = numeric_limits<Real>::tinyness_before;
    static constexpr float_round_style round_style = numeric_limits<Real>::round_style;
  };
}
