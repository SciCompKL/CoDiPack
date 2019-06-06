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

#include "dummyValue.hpp"
#include "vectorStorage.hpp"
#include "../../configure.h"

/**
 * @brief Global namespace for CoDiPack - Code Differentiation Package
 */
namespace codi {

  struct DummyHessian {
      CODI_INLINE DummyValue operator()(const size_t i, const size_t j, const size_t k) {
        CODI_UNUSED(i);
        CODI_UNUSED(j);
        CODI_UNUSED(k);

        return DummyValue();
      }
  };

  template <typename Vec>
  struct Hessian {

      VectorStorage<Vec> values;
      using T = typename VectorStorage<Vec>::Element;

      size_t m;
      size_t n;

      explicit Hessian(size_t m, size_t n) : values(n * n * m), m(m), n(n) {}

      CODI_INLINE T operator()(const size_t i, const size_t j, const size_t k) const {
        return values.data()[computeIndex(i,j,k)];
      }

      CODI_INLINE T& operator()(const size_t i, const size_t j, const size_t k) {
        return values.data()[computeIndex(i,j,k)];
      }

      template<typename T>
      CODI_INLINE void set(const size_t i, const size_t j, const size_t k, const T& v) {
        values.data()[computeIndex(i,j,k)] = v;
      }

    private:

      CODI_INLINE size_t computeIndex(const size_t i, const size_t j, const size_t k) const {
        return k * n * m + i * n + j;
      }
  };

}
