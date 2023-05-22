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
#include "../../misc/macros.hpp"
#include "../../traits/realTraits.hpp"
#include "delayAccessor.hpp"
#include "dummy.hpp"
#include "jacobianInterface.hpp"
#include "staticDummy.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Default implementation of the Jacobian interface.
   *
   * Data is stored in a row-major format.
   *
   * @tparam T_T  The data type in the Jacobian.
   * @tparam T_Store  Storage allocator. Should implement the standard vector interface.
   */
  template<typename T_T, typename T_Store = std::vector<T_T>>
  struct Jacobian : public JacobianInterface<T_T> {
    public:

      using T = CODI_DD(T_T, double);                       ///< See Jacobian.
      using Store = CODI_DD(T_Store, std::vector<double>);  ///< See Jacobian.

    protected:
      Store values;  ///< Array store for the data.

      size_t m;  ///< Number of rows (output variables).
      size_t n;  ///< Number of columns (input variables).

    public:

      /// m = rows (output variables), n = columns (input variables)
      explicit Jacobian(size_t const m, size_t const n) : values(constructVector<Store>(n * m)), m(m), n(n) {}

      /// \copydoc codi::JacobianInterface::getM()
      CODI_INLINE size_t getM() const {
        return m;
      }

      /// \copydoc codi::JacobianInterface::getN()
      CODI_INLINE size_t getN() const {
        return n;
      }

      /// \copydoc codi::JacobianInterface::operator()(size_t const, size_t const) const
      CODI_INLINE T operator()(size_t const i, size_t const j) const {
        return values[computeIndex(i, j)];
      }

      /// \copydoc codi::JacobianInterface::operator()(size_t const, size_t const)
      CODI_INLINE T& operator()(size_t const i, size_t const j) {
        return values[computeIndex(i, j)];
      }

      /// \copydoc codi::JacobianInterface::resize()
      ///
      /// Implementation: Old values are not overwritten.
      CODI_INLINE void resize(size_t const m, size_t const n) {
        this->m = m;
        this->n = n;

        values.resize(m * n);
      }

      /// \copydoc codi::JacobianInterface::size()
      CODI_INLINE size_t size() const {
        return m * n;
      }

    protected:

      /// Convert row and column to an index in the storage array.
      CODI_INLINE size_t computeIndex(size_t const i, size_t const j) const {
        return i * n + j;
      }
  };

  /**
   * @brief Adds counting of nonzero entries.
   *
   * The user has to manually reset the count.
   *
   * @tparam T_T  The data type in the Jacobian.
   * @tparam T_Store  Storage allocator. Should implement the standard vector interface.
   */
  template<typename T_T, typename T_Store = std::vector<T_T>>
  struct JacobianCountNonZerosRow : public Jacobian<T_T, T_Store> {
    public:

      using Base = Jacobian<T_T, T_Store>;  ///< Base class abbreviation.

      using T = typename Base::T;          ///< See Jacobian.
      using Store = typename Base::Store;  ///< See Jacobian.

      using DelayAcc = JacobianDelayAccessor<JacobianCountNonZerosRow>;  ///< Delayed accessor for reference access.

    private:

      std::vector<int> nonZerosRowVector;  ///< Count nonzero entries in each row.

    public:

      /// \copydoc codi::Jacobian::Jacobian
      explicit JacobianCountNonZerosRow(size_t m, size_t n) : Base(m, n), nonZerosRowVector(m) {}

      /// \copydoc codi::JacobianInterface::operator()(size_t const  i, size_t const j) const
      CODI_INLINE T operator()(size_t const i, size_t const j) const {
        return Base::operator()(i, j);
      }

      /// \copydoc codi::JacobianInterface::operator()(size_t const i, size_t const j)
      /// Implementation: Returns an object for the delayed access. This object calls then the access logic which
      /// updates the number of nonzero elements.
      CODI_INLINE DelayAcc operator()(size_t const i, size_t const j) {
        return DelayAcc(i, j, *this);
      }

      /// \copydoc codi::JacobianInterface::resize
      CODI_INLINE void resize(size_t const m, size_t const n) {
        Base::resize(m, n);
        nonZerosRowVector.resize(m);
      }

      /// Reference to the number of nonzero entries for the specified row.
      CODI_INLINE int& nonZerosRow(size_t const i) {
        return nonZerosRowVector[i];
      }

      /// Checks if the element is non zero and updates the counter.
      CODI_INLINE void setLogic(size_t const i, size_t const j, T const& v) {
        if (T() != v) {
          nonZerosRowVector[i] += 1;
        }
        Base::operator()(i, j) = v;
      }
  };

  /**
   * @brief Wrapper for JacboianInterfaces that requires a passive value conversion.
   *
   * @tparam T_Nested  The nested Jacobian, that has a passive value storage.
   */
  template<typename T_Nested>
  struct JacobianConvertWrapper {
    public:
      using Nested = CODI_DD(T_Nested, JacobianInterface<double>);     ///< See JacobianConvertWrapper.
      using T = typename Nested::T;                                    ///< See JacobianInterface.
      using DelayAcc = JacobianDelayAccessor<JacobianConvertWrapper>;  ///< Return type for reference access.

    private:
      Nested& nested;

    public:

      /// Constructor
      explicit JacobianConvertWrapper(Nested& nested) : nested(nested) {}

      /// \copydoc JacobianInterface::operator()(size_t const  i, size_t const j) const
      CODI_INLINE T operator()(size_t const i, size_t const j) const {
        return nested(i, j);
      }

      /// \copydoc JacobianInterface::operator()(size_t const  i, size_t const j)
      /// Returns a JacobianDelayAccessor for the tracking of the nonzero entries.
      CODI_INLINE DelayAcc operator()(size_t const i, size_t const j) {
        return DelayAcc(i, j, *this);
      }

      /// Called by the JacobianDelayAccessor when a value is set.
      template<typename SetT>
      CODI_INLINE void setLogic(size_t const i, size_t const j, SetT const& v) {
        nested(i, j) = RealTraits::getPassiveValue<SetT>(v);
      }
  };

  /// Dummy Jacobian. Has size zero and no logic in any call.
  struct DummyJacobian : public JacobianInterface<DummyValue> {
    public:

      /// \copydoc JacobianInterface::getM()
      size_t getM() const {
        return 0;
      }

      /// \copydoc JacobianInterface::getN()
      size_t getN() const {
        return 0;
      }

      /// \copydoc JacobianInterface::operator()(size_t const  i, size_t const j) const
      CODI_INLINE DummyValue operator()(size_t const i, size_t const j) const {
        CODI_UNUSED(i, j);

        return DummyValue();
      }

      /// \copydoc JacobianInterface::operator()(size_t const  i, size_t const j)
      CODI_INLINE DummyValue& operator()(size_t const i, size_t const j) {
        CODI_UNUSED(i, j);

        return StaticDummy<DummyValue>::dummy;
      }

      /// \copydoc JacobianInterface::resize()
      void resize(size_t const m, size_t const n) {
        CODI_UNUSED(m, n);
      }

      /// \copydoc JacobianInterface::size()
      size_t size() const {
        return 0;
      }
  };
}
