#pragma once

#include <vector>

#include "../../config.h"
#include "dummy.hpp"
#include "staticDummy.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  template<typename _T>
  struct HessianInterface {
    public:

      using T = CODI_DECLARE_DEFAULT(_T, double);

      virtual ~HessianInterface() {}

      size_t getM() const;
      size_t getN() const;

      T operator()(size_t const i, size_t const j, size_t const k) const;
      T& operator()(size_t const i, size_t const j, size_t const k);

      void resize(size_t const m, size_t const n);
      size_t size() const;
  };

  template <typename _T, typename _Store = std::vector<_T>>
  struct Hessian : public HessianInterface<_T> {
    public:

      using T = CODI_DECLARE_DEFAULT(_T, double);
      using Store = CODI_DECLARE_DEFAULT(_Store, std::vector<double>);

    protected:
      Store values;

      size_t m;
      size_t n;

    public:
      explicit Hessian(size_t m, size_t n) : values(n * n * m), m(m), n(n) {}

      CODI_INLINE size_t getM() const {return m;}
      CODI_INLINE size_t getN() const {return n;}

      CODI_INLINE T operator()(size_t const i, size_t const j, size_t const k) const {
        return values.data()[computeIndex(i,j,k)];
      }

      CODI_INLINE T& operator()(size_t const i, size_t const j, size_t const k) {
        return values.data()[computeIndex(i,j,k)];
      }

      CODI_INLINE void resize(size_t const m, size_t const n) {
        this->m = m;
        this->n = n;

        values.resize(m * n * n);
      }

      CODI_INLINE size_t size() const {return m * n * n;}

    private:
      CODI_INLINE size_t computeIndex(size_t const i, size_t const j, size_t const k) const {
        return k * n * m + i * n + j;
      }
  };

  struct DummyHessian : public HessianInterface<DummyValue> {
    public:

      size_t getM() const {return 0;}
      size_t getN() const {return 0;}

      CODI_INLINE DummyValue operator()(size_t const i, size_t const j, size_t const k) const {
        CODI_UNUSED(i, j, k);

        return DummyValue();
      }

      CODI_INLINE DummyValue& operator()(size_t const i, size_t const j, size_t const k) {
        CODI_UNUSED(i, j, k);

        return StaticDummy<DummyValue>::dummy;
      }

      void resize(size_t const m, size_t const n) {CODI_UNUSED(m, n);}
      size_t size() const {return 0;}
  };
}
