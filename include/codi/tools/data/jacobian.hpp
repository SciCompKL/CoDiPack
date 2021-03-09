#pragma once

#include <vector>

#include "../../aux/macros.hpp"
#include "../../config.h"
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
   * @tparam _T  The data type in the Jacobian.
   * @tparam _Store  Storage allocator. Should implement the standard vector interface.
   */
  template<typename _T, typename _Store = std::vector<_T>>
  struct Jacobian : public JacobianInterface<_T> {
    public:

      using T = CODI_DD(_T, double);                       ///< See Jacobian
      using Store = CODI_DD(_Store, std::vector<double>);  ///< See Jacobian

    protected:
      Store values;  ///< Array store for the data.

      size_t m;  ///< Number of row (output variables)
      size_t n;  ///< Number of columns (input variables);

    public:

      /// m = rows ( output variables), n = columns (input variables)
      explicit Jacobian(size_t const m, size_t const n) : values(n * m), m(m), n(n) {}

      CODI_INLINE size_t getM() const {
        return m;
      }  ///< \copydoc codi::JacobianInterface::getM()
      CODI_INLINE size_t getN() const {
        return n;
      }  ///< \copydoc codi::JacobianInterface::getN()

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

      /// Compute index into the storage array.
      CODI_INLINE size_t computeIndex(size_t const i, size_t const j) const {
        return i * n + j;
      }
  };

  template<typename _T, typename _Store = std::vector<_T>>
  struct JacobianCountNonZerosRow : public Jacobian<_T, _Store> {
    public:

      using Base = Jacobian<_T, _Store>;  ///< Base class abbreviation

      using T = typename Base::T;          ///< See Jacobian
      using Store = typename Base::Store;  ///< See Jacobian

      using DelayAcc = DelayAccessor<JacobianCountNonZerosRow>;

    private:

      std::vector<int> nonZerosRowVector;  ///< Count nonzero entries in each row.

    public:

      /// \copydoc Jacobian::Jacobian
      explicit JacobianCountNonZerosRow(size_t m, size_t n) : Base(m, n), nonZerosRowVector(m) {}

      /// \copydoc JacobianInterface::operator()(size_t const  i, size_t const j) const
      CODI_INLINE T operator()(size_t const i, size_t const j) const {
        return Base::operator()(i, j);
      }

      /// \copydoc JacobianInterface::operator()(size_t const i, size_t const j)
      /// Implementation: Returns an object for the delayed access. This object calls then the access logic which
      /// updates the number of non zero elements.
      CODI_INLINE DelayAcc operator()(size_t const i, size_t const j) {
        return DelayAcc(i, j, *this);
      }

      /// \copydoc JacobianInterface::resize
      CODI_INLINE void resize(size_t const m, size_t const n) {
        Base::resize(m, n);
        nonZerosRowVector.resize(m);
      }

      /// \copydoc JacobianInterface::reshape
      CODI_INLINE void reshape(size_t const m, size_t const n) {
        Base::reshape(m, n);
        if (nonZerosRowVector.size() < m) {
          nonZerosRowVector.resize(m);
        }
      }

      /// Reference to the number of non zero entries for the specified row.
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

  template<typename _Nested>
  struct JacobianConvertWrapper {
    public:
      using Nested = CODI_DECLARE_DEFAULT(_Nested, JacobianInterface<double>);
      using T = typename Nested::T;
      using DelayAcc = DelayAccessor<JacobianConvertWrapper>;

    private:
      Nested& nested;

    public:
      explicit JacobianConvertWrapper(Nested& nested) : nested(nested) {}

      CODI_INLINE T operator()(size_t const i, size_t const j) const {
        return nested(i, j);
      }

      CODI_INLINE DelayAcc operator()(size_t const i, size_t const j) {
        return DelayAcc(i, j, *this);
      }

      template<typename SetT>
      CODI_INLINE void setLogic(size_t const i, size_t const j, SetT const& v) {
        nested(i, j) = RealTraits::getPassiveValue<SetT>(v);
      }
  };

  struct DummyJacobian : public JacobianInterface<DummyValue> {
    public:

      size_t getM() const {
        return 0;
      }
      size_t getN() const {
        return 0;
      }

      CODI_INLINE DummyValue operator()(size_t const i, size_t const j) const {
        CODI_UNUSED(i, j);

        return DummyValue();
      }

      CODI_INLINE DummyValue& operator()(size_t const i, size_t const j) {
        CODI_UNUSED(i, j);

        return StaticDummy<DummyValue>::dummy;
      }

      void resize(size_t const m, size_t const n) {
        CODI_UNUSED(m, n);
      }
      size_t size() const {
        return 0;
      }
  };
}
