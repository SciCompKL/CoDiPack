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

#include <vector>
#include <array>

#include "../../configure.h"
#include "../../macros.h"

/**
 * @brief Global namespace for CoDiPack - Code Differentiation Package
 */
namespace codi {

  template <typename Vec>
  struct VectorStorage {
      using VecType = Vec;
      using Element = char;

      explicit VectorStorage(size_t size);

      CODI_INLINE Element const* data() const;

      CODI_INLINE Element* data();

      CODI_INLINE Element& operator[](size_t i);

      CODI_INLINE size_t size();

      CODI_NO_INLINE void resize(size_t const size);
  };

  template <typename T, typename Allocator>
  struct VectorStorage<std::vector<T, Allocator>> {

      using VecType = std::vector<T, Allocator>;
      using Element = T;

      VecType vec;

      explicit VectorStorage(size_t size) : vec(size) {}

      CODI_INLINE T const* data() const {
        return vec.data();
      }

      CODI_INLINE T* data() {
        return vec.data();
      }

      CODI_INLINE T& operator[](size_t i) {
        return vec[i];
      }

      CODI_INLINE size_t size() {
        return vec.size();
      }

      CODI_NO_INLINE void resize(size_t const size) {
        vec.resize(size);
      }
  };

  template <typename T, size_t N>
  struct VectorStorage<std::array<T, N>> {

      using VecType = std::array<T, N>;
      using Element = T;


      VecType vec;

      explicit VectorStorage(size_t size) : vec() {CODI_UNUSED(size);}

      CODI_INLINE T const* data() const {
        return vec.data();
      }

      CODI_INLINE T* data() {
        return vec.data();
      }

      CODI_INLINE T& operator[](size_t i) {
        return vec[i];
      }

      CODI_INLINE size_t size() {
        return N;
      }

      CODI_INLINE void resize(size_t const size) {
        CODI_UNUSED(size);

        CODI_EXCEPTION("Can not resize std::array.");
      }
  };
}
