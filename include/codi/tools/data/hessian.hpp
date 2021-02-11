#pragma once

#include "../../config.h"

#include <vector>

/** \copydoc codi::Namespace */
namespace codi {

  template<typename _T>
  struct HessianInterface {
    public:

      using T = CODI_DECLARE_DEFAULT(_T, double);

      virtual ~HessianInterface() {}

      virtual T operator()(size_t const i, size_t const j, size_t const k) const = 0;
      virtual T& operator()(size_t const i, size_t const j, size_t const k) = 0;
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

      template<typename T>
      CODI_INLINE void set(size_t const i, size_t const j, size_t const k, T& const v) {
        values.data()[computeIndex(i,j,k)] = v;
      }

    private:

      CODI_INLINE size_t computeIndex(size_t const i, size_t const j, size_t const k) const {
        return k * n * m + i * n + j;
      }
  };

}
