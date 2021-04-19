#pragma once

#include <vector>

#include "../../aux/constructVector.hpp"
#include "../../config.h"
#include "dummy.hpp"
#include "hessianInterface.hpp"
#include "staticDummy.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Default implementation of the Hessian interface.
   *
   * Running index speed: j(fastest), i, k(slowest)
   *
   * Data is stored in an array of row-major matrices.
   *
   * @tparam _T  Data type in the Hessian.
   * @tparam _Store  Storage allocator. Should implement the standard vector interface.
   */
  template<typename _T, typename _Store = std::vector<_T>>
  struct Hessian : public HessianInterface<_T> {
    public:

      using T = CODI_DECLARE_DEFAULT(_T, double);                       ///< See Hessian
      using Store = CODI_DECLARE_DEFAULT(_Store, std::vector<double>);  ///< See Hessian

    protected:
      Store values;  ///< Storage for the data

      size_t m;  ///< Number of function outputs
      size_t n;  ///< Number of function inputs

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

  /// Dummy Hessian. Has size zero and no logic in all calls.
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
