#pragma once

#include <vector>

#include "../../aux/macros.hpp"
#include "../../config.h"
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
  template<typename _T, typename _Store = std::vector<_T>>
  struct Jacobian : public JacobianInterface<_T> {
    public:

      using T = CODI_DD(_T, double);  ///< See Jacobian
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

      /// \copydoc codi::JacobianInterface::operator()( size_t const, size_t const) const
      CODI_INLINE T operator()(size_t const i, size_t const j) const {
        return values[computeIndex(i, j)];
      }

      /// \copydoc codi::JacobianInterface::operator()( size_t const, size_t const)
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

      ///< Compute index into the storage array.
      CODI_INLINE size_t computeIndex(const size_t i, const size_t j) const {
        return i * n + j;
      }
  };
}
