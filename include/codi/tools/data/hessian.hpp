/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2023 Chair for Scientific Computing (SciComp), University of Kaiserslautern-Landau
 * Homepage: http://www.scicomp.uni-kl.de
 * Contact:  Prof. Nicolas R. Gauger (codi@scicomp.uni-kl.de)
 *
 * Lead developers: Max Sagebaum, Johannes Blühdorn (SciComp, University of Kaiserslautern-Landau)
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
 * For other licensing options please contact us.
 *
 * Authors:
 *  - SciComp, University of Kaiserslautern-Landau:
 *    - Max Sagebaum
 *    - Johannes Blühdorn
 *    - Former members:
 *      - Tim Albring
 */
#pragma once

#include <vector>

#include "../../config.h"
#include "../../misc/constructVector.hpp"
#include "dummy.hpp"
#include "hessianInterface.hpp"
#include "staticDummy.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Default implementation of the Hessian interface.
   *
   * Running index speed: j (fastest), i, k (slowest).
   *
   * Data is stored in an array of row-major matrices.
   *
   * @tparam T_T  Data type in the Hessian.
   * @tparam T_Store  Storage allocator. Should implement the standard vector interface.
   */
  template<typename T_T, typename T_Store = std::vector<T_T>>
  struct Hessian : public HessianInterface<T_T> {
    public:

      using T = CODI_DD(T_T, double);                       ///< See Hessian.
      using Store = CODI_DD(T_Store, std::vector<double>);  ///< See Hessian.

    protected:
      Store values;  ///< Storage for the data.

      size_t m;  ///< Number of function outputs.
      size_t n;  ///< Number of function inputs.

    public:

      /// Constructor
      explicit Hessian(size_t m, size_t n) : values(constructVector<Store>(n * n * m)), m(m), n(n) {}

      /// \copydoc HessianInterface::getM()
      CODI_INLINE size_t getM() const {
        return m;
      }

      /// \copydoc HessianInterface::getN()
      CODI_INLINE size_t getN() const {
        return n;
      }

      /// \copydoc HessianInterface::operator()(size_t const i, size_t const j, size_t const k) const
      CODI_INLINE T operator()(size_t const i, size_t const j, size_t const k) const {
        return values.data()[computeIndex(i, j, k)];
      }

      /// \copydoc HessianInterface::operator()(size_t const i, size_t const j, size_t const k)
      CODI_INLINE T& operator()(size_t const i, size_t const j, size_t const k) {
        return values.data()[computeIndex(i, j, k)];
      }

      /// \copydoc HessianInterface::resize()
      CODI_INLINE void resize(size_t const m, size_t const n) {
        this->m = m;
        this->n = n;

        values.resize(m * n * n);
      }

      /// \copydoc HessianInterface::size()
      CODI_INLINE size_t size() const {
        return m * n * n;
      }

    private:
      CODI_INLINE size_t computeIndex(size_t const i, size_t const j, size_t const k) const {
        return k * n * m + i * n + j;
      }
  };

  /// Dummy Hessian. Has size zero and no logic in any call.
  struct DummyHessian : public HessianInterface<DummyValue> {
    public:

      /// \copydoc HessianInterface::getM()
      size_t getM() const {
        return 0;
      }

      /// \copydoc HessianInterface::getN()
      size_t getN() const {
        return 0;
      }

      /// \copydoc HessianInterface::operator()(size_t const i, size_t const j, size_t const k) const
      CODI_INLINE DummyValue operator()(size_t const i, size_t const j, size_t const k) const {
        CODI_UNUSED(i, j, k);

        return DummyValue();
      }

      /// \copydoc HessianInterface::operator()(size_t const i, size_t const j, size_t const k)
      CODI_INLINE DummyValue& operator()(size_t const i, size_t const j, size_t const k) {
        CODI_UNUSED(i, j, k);

        return StaticDummy<DummyValue>::dummy;
      }

      /// \copydoc HessianInterface::resize()
      void resize(size_t const m, size_t const n) {
        CODI_UNUSED(m, n);
      }

      /// \copydoc HessianInterface::size()
      size_t size() const {
        return 0;
      }
  };
}
