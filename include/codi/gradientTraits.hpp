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

#include "configure.h"
#include "macros.h"
#include "tools/direction.hpp"

/**
 * @brief Global namespace for CoDiPack - Code Differentiation Package
 */
namespace codi {

  template<typename T>
  struct GradientValueTraits {
      using Data = T;

      static CODI_INLINE constexpr size_t getVectorSize() {
        return 1;
      }

      static CODI_INLINE T& at(T& v, const size_t pos) {
        CODI_UNUSED(pos);

        return v;
      }

      static CODI_INLINE const T& at(const T& v, const size_t pos) {
        CODI_UNUSED(pos);

        return v;
      }
  };

  template<typename T, size_t n>
  struct GradientValueTraits<Direction<T, n>> {

      using Data = T;
      static CODI_INLINE constexpr size_t getVectorSize() {
        return n;
      }

      static CODI_INLINE T& at(Direction<T,n>& v, const size_t pos) {
        return v[pos];
      }

      static CODI_INLINE const T& at(const Direction<T,n>& v, const size_t pos) {
        return v[pos];
      }
  };
}
