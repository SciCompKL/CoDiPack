#pragma once

#include <vector>

#include "../../aux/macros.hpp"
#include "../../config.h"
#include "jacobianInterface.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  template <typename _T, typename _Store = std::vector<_T>>
  struct Jacobian : public JacobianInterface<_T>{
    public:

      using T = CODI_DECLARE_DEFAULT(_T, double);
      using Store = CODI_DECLARE_DEFAULT(_Store, std::vector<double>);

    private:
      Store values;

      size_t m;
      size_t n;

    public:

      explicit Jacobian(size_t const m, size_t const n) : values(n * m), m(m), n(n) {}

      CODI_INLINE size_t getM() const {return m;}
      CODI_INLINE size_t getN() const {return n;}

      CODI_INLINE T operator()(size_t const i, size_t const j) const {
        return values[computeIndex(i,j)];
      }

      CODI_INLINE T& operator()(size_t const i, size_t const j) {
        return values[computeIndex(i,j)];
      }

      CODI_INLINE void resize(size_t const m, size_t const n) {
        this->m = m;
        this->n = n;

        values.resize(m*n);
      }

      CODI_INLINE size_t size() const {return m * n;}

    protected:

      CODI_INLINE size_t computeIndex(const size_t i, const size_t j) const {
        return i * n + j;
      }
  };
}
