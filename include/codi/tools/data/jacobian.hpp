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

/**
 * @brief Global namespace for CoDiPack - Code Differentiation Package
 */
namespace codi {

  struct DummyJacobian {
      CODI_INLINE DummyValue operator()(const size_t i, const size_t j) {
        CODI_UNUSED(i);
        CODI_UNUSED(j);

        return DummyValue();
      }
  };

  template <typename Vec>
  struct Jacobian {

      VectorStorage<Vec> values;
      using T = typename VectorStorage<Vec>::Element;

      size_t m;
      size_t n;

      explicit Jacobian(size_t m, size_t n) : values(n * m), m(m), n(n) {}

      CODI_INLINE T operator()(const size_t i, const size_t j) const {
        return values.data()[computeIndex(i,j)];
      }

      CODI_INLINE T& operator()(const size_t i, const size_t j) {
        return values.data()[computeIndex(i,j)];
      }

      CODI_INLINE void resize(size_t const m, size_t const n) {
        this->m = m;
        this->n = n;

        values.resize(m*n);
      }

    protected:

      CODI_INLINE size_t computeIndex(const size_t i, const size_t j) const {
        return i * n + j;
      }
  };

  template <typename Vec>
  struct JacobianCountNonZerosRow {

      using T = typename VectorStorage<Vec>::Element;

      struct ValueAccessor {
          size_t i;
          size_t j;

          JacobianCountNonZerosRow& data;

          ValueAccessor(size_t const i, size_t const j, JacobianCountNonZerosRow& data) : i(i), j(j), data(data) {}

          ValueAccessor& operator =(JacobianCountNonZerosRow::T const& v) {
            data.set(i,j, v);

            return *this;
          }

          operator T() const {
            return data.get(i, j);
          }
      };

      VectorStorage<Vec> values;
      std::vector<int> nonZerosRow;

      size_t m;
      size_t n;

      explicit JacobianCountNonZerosRow(size_t m, size_t n) : values(n * m), nonZerosRow(m), m(m), n(n) {}

      CODI_INLINE T operator()(const size_t i, const size_t j) const {
        return get(i,j);
      }

      CODI_INLINE ValueAccessor operator()(const size_t i, const size_t j) {
        return ValueAccessor(i,j,*this);
      }

      CODI_INLINE T get(const size_t i, const size_t j) const {
        return values.data()[computeIndex(i,j)];
      }

      CODI_INLINE void set(const size_t i, const size_t j, T const& v) {
        if(T() != v) {
          nonZerosRow[i] += 1;
          values.data()[computeIndex(i,j)] = v;
        }
      }

      CODI_INLINE void resize(size_t const m, size_t const n) {
        this->m = m;
        this->n = n;

        values.resize(m*n);
        nonZerosRow.resize(m);
      }

      CODI_INLINE size_t size() {
        return m * n;
      }

    protected:

      CODI_INLINE size_t computeIndex(const size_t i, const size_t j) const {
        return i * n + j;
      }
  };

}
