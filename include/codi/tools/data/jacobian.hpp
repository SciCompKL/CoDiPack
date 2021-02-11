#pragma once

#include <vector>

#include "../../aux/macros.hpp"
#include "../../config.h"
#include "../../traits/realTraits.hpp"
#include "jacobianInterface.hpp"

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
  template <typename _T, typename _Store = std::vector<_T>>
  struct Jacobian : public JacobianInterface<_T>{
    public:

      using T = CODI_DECLARE_DEFAULT(_T, double); ///< See Jacobian
      using Store = CODI_DECLARE_DEFAULT(_Store, std::vector<double>); ///< See Jacobian

    protected:
      Store values; ///< Array store for the data.

      size_t m; ///< Number of row (output variables)
      size_t n; ///< Number of columns (input variables);

    public:

      /// m = rows ( output variables), n = columns (input variables)
      explicit Jacobian(size_t const m, size_t const n) : values(n * m), m(m), n(n) {}

      CODI_INLINE size_t getM() const {return m;} ///< \copydoc codi::JacobianInterface::getM()
      CODI_INLINE size_t getN() const {return n;} ///< \copydoc codi::JacobianInterface::getN()

      /// \copydoc codi::JacobianInterface::operator()(size_t const, size_t const) const
      CODI_INLINE T operator()(size_t const i, size_t const j) const {
        return values[computeIndex(i,j)];
      }

      /// \copydoc codi::JacobianInterface::operator()(size_t const, size_t const)
      CODI_INLINE T& operator()(size_t const i, size_t const j) {
        return values[computeIndex(i,j)];
      }

      /// \copydoc codi::JacobianInterface::resize()
      ///
      /// Implementation: Old values are not overwritten.
      CODI_INLINE void resize(size_t const m, size_t const n) {
        this->m = m;
        this->n = n;

        values.resize(m*n);
      }

      /// \copydoc codi::JacobianInterface::size()
      CODI_INLINE size_t size() const {return m * n;}

    protected:

      ///< Compute index into the storage array.
      CODI_INLINE size_t computeIndex(size_t const i, size_t const j) const {
        return i * n + j;
      }
  };


  template <typename Impl>
  struct DelayAccessor {
    private:
      size_t i;
      size_t j;

      Impl& data;

    public:

      DelayAccessor(size_t const i, size_t const j, Impl& data) : i(i), j(j), data(data) {}

      template<typename T>
      DelayAccessor& operator =(T const& v) {
        data.setLogic(i,j, v);

        return *this;
      }

      operator typename Impl::T() const {
        return const_cast<Impl const&>(data).operator()(i,j);
      }
  };


  template <typename _T, typename _Store = std::vector<_T>>
  struct JacobianCountNonZerosRow : public Jacobian<_T, _Store> {
    public:

      std::vector<int> nonZerosRow;

    public:
      using T = typename Jacobian<_T, _Store>::T; ///< See Jacobian
      using Store = typename Jacobian<_T, _Store>::Store; ///< See Jacobian

      using DelayAcc = DelayAccessor<JacobianCountNonZerosRow>;

      /** \copydoc Jacobian::Jacobian */
      explicit JacobianCountNonZerosRow(size_t m, size_t n) : Jacobian<T, Store>(m, n), nonZerosRow(m) {}


      /** \copydoc Jacobian::operator()(size_t const  i, size_t const j) const */
      CODI_INLINE T operator()(size_t const i, size_t const j) const {
        return Jacobian<T, Store>::operator()(i, j);
      }

      /**
       * @brief Returns an object for the delayed access. This object calls then the access logic which updates the
       * number of non zero elements.
       *
       * \copydetails Jacobian::operator()(size_t const i, size_t const j)
       */
      CODI_INLINE DelayAcc operator()(size_t const i, size_t const j) {
        return DelayAcc(i, j, *this);
      }

      /** \copydoc Jacobian::resize */
      CODI_INLINE void resize(size_t const m, size_t const n) {

        Jacobian<T, Store>::resize(m, n);
        nonZerosRow.resize(m);
      }

      /** \copydoc Jacobian::reshape */
      CODI_INLINE void reshape(size_t const m, size_t const n) {
        Jacobian<T, Store>::reshape(m, n);
        if (nonZerosRow.size() < m) {
          nonZerosRow.resize(m);
        }
      }
      /**
       * @brief Reference to the number of non zero entries for the specified row.
       *
       * @param[in] i  Output value position. Range: [0, m)
       *
       * @return Reference to the non zero count for the specified row.
       */
      CODI_INLINE int& nonZeroRow(size_t const i) {
        return nonZerosRow[i];
      }

      /**
       * @brief Checks if the element is non zero and updates the counter.
       *
       * @param[in] i  Default
       * @param[in] j  Default
       * @param[in] v  Default
       */
      CODI_INLINE void setLogic(size_t const i, size_t const j, T const& v) {
        if (T() != v) {
          nonZerosRow[i] += 1;
        }
        Jacobian<T, Store>::operator()(i, j) = v;
      }
  };

  template <typename _Nested>
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

}
