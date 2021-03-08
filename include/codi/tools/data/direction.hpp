#pragma once

#include <initializer_list>

#include "../../aux/macros.hpp"
#include "../../config.h"
#include "../../traits/gradientTraits.hpp"
#include "../../traits/realTraits.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Fixed size vector mode implementation.
   *
   * Can be used as the gradient template argument in active CoDiPack types.
   *
   * @tparam _Real  Value of the gradients
   * @tparam _dim  Dimension of the vector mode.
   */
  template<typename _Real, size_t _dim>
  struct Direction {
    public:

      using Real = CODI_DD(_Real, double);  ///< See Direction

      static size_t constexpr dim = _dim;  ///< See Direction

    private:
      Real vector[dim];

    public:

      /// Constructor
      CODI_INLINE Direction() : vector() {}

      /// Constructor
      CODI_INLINE Direction(Real const& s) : vector() {
        for (size_t i = 0; i < dim; ++i) {
          vector[i] = s;
        }
      }

      /// Constructor
      CODI_INLINE Direction(Direction const& v) : vector() {
        for (size_t i = 0; i < dim; ++i) {
          vector[i] = v.vector[i];
        }
      }

      /// Constructor
      CODI_INLINE Direction(std::initializer_list<Real> l) : vector() {
        size_t size = std::min(dim, l.size());
        Real const* array = l.begin();  // this is possible because the standard requires an array storage
        for (size_t i = 0; i < size; ++i) {
          vector[i] = array[i];
        }
      }

      /// Per reference element access.
      CODI_INLINE Real& operator[](size_t const& i) {
        return vector[i];
      }

      /// Per value element access.
      CODI_INLINE Real const& operator[](size_t const& i) const {
        return vector[i];
      }

      /// Assignment operator
      CODI_INLINE Direction& operator=(Direction const& v) {
        for (size_t i = 0; i < dim; ++i) {
          this->vector[i] = v.vector[i];
        }

        return *this;
      }

      /// Update operator
      CODI_INLINE Direction& operator+=(Direction const& v) {
        for (size_t i = 0; i < dim; ++i) {
          this->vector[i] += v.vector[i];
        }

        return *this;
      }

      /// Update operator
      CODI_INLINE Direction& operator-=(Direction const& v) {
        for (size_t i = 0; i < dim; ++i) {
          this->vector[i] -= v.vector[i];
        }

        return *this;
      }
  };

  template<typename Real, size_t dim>
  size_t constexpr Direction<Real, dim>::dim;

  /// Multiplication with scalar
  template<typename Real, size_t dim>
  CODI_INLINE Direction<Real, dim> operator*(Real const& s, Direction<Real, dim> const& v) {
    Direction<Real, dim> r;
    for (size_t i = 0; i < dim; ++i) {
      r[i] = s * v[i];
    }

    return r;
  }

  /// Multiplication with passive scalar
  template<typename Real, size_t dim, typename = RealTraits::EnableIfNotPassiveReal<Real>>
  CODI_INLINE Direction<Real, dim> operator*(RealTraits::PassiveReal<Real> const& s, Direction<Real, dim> const& v) {
    Direction<Real, dim> r;
    for (size_t i = 0; i < dim; ++i) {
      r[i] = s * v[i];
    }

    return r;
  }

  /// Multiplication with scalar
  template<typename Real, size_t dim>
  CODI_INLINE Direction<Real, dim> operator*(Direction<Real, dim> const& v, Real const& s) {
    return s * v;
  }

  /// Multiplication with passive scalar
  template<typename Real, size_t dim, typename = RealTraits::EnableIfNotPassiveReal<Real>>
  CODI_INLINE Direction<Real, dim> operator*(Direction<Real, dim> const& v, RealTraits::PassiveReal<Real> const& s) {
    return s * v;
  }

  /// Division with scalar
  template<typename Real, size_t dim>
  CODI_INLINE Direction<Real, dim> operator/(Direction<Real, dim> const& v, Real const& s) {
    Direction<Real, dim> r;
    for (size_t i = 0; i < dim; ++i) {
      r[i] = v[i] / s;
    }

    return r;
  }

  /// Division with passive scalar
  template<typename Real, size_t dim, typename = RealTraits::EnableIfNotPassiveReal<Real>>
  CODI_INLINE Direction<Real, dim> operator/(Direction<Real, dim> const& v, RealTraits::PassiveReal<Real> const& s) {
    Direction<Real, dim> r;
    for (size_t i = 0; i < dim; ++i) {
      r[i] = v[i] / s;
    }

    return r;
  }

  /// Summation of two vectors
  template<typename Real, size_t dim>
  CODI_INLINE Direction<Real, dim> operator+(Direction<Real, dim> const& v1, Direction<Real, dim> const& v2) {
    Direction<Real, dim> r;
    for (size_t i = 0; i < dim; ++i) {
      r[i] = v1[i] + v2[i];
    }

    return r;
  }

  /// Subtraction of two vectors
  template<typename Real, size_t dim>
  CODI_INLINE Direction<Real, dim> operator-(Direction<Real, dim> const& v1, Direction<Real, dim> const& v2) {
    Direction<Real, dim> r;
    for (size_t i = 0; i < dim; ++i) {
      r[i] = v1[i] - v2[i];
    }

    return r;
  }

  /// Negation
  template<typename Real, size_t dim>
  CODI_INLINE Direction<Real, dim> operator-(Direction<Real, dim> const& v) {
    Direction<Real, dim> r;
    for (size_t i = 0; i < dim; ++i) {
      r[i] = -v[i];
    }

    return r;
  }

  /// True if one element differs.
  template<typename A, typename Real, size_t dim>
  CODI_INLINE bool operator!=(A const& s, Direction<Real, dim> const& v) {
    for (size_t i = 0; i < dim; ++i) {
      if (s != v[i]) {
        return true;
      }
    }

    return false;
  }

  /// True if one element differs.
  template<typename A, typename Real, size_t dim>
  CODI_INLINE bool operator!=(Direction<Real, dim> const& v, A const& s) {
    return s != v;
  }

  /// True if all elements are the same.
  template<typename A, typename Real, size_t dim>
  CODI_INLINE bool operator==(A const& s, Direction<Real, dim> const& v) {
    for (size_t i = 0; i < dim; ++i) {
      if (s != v[i]) {
        return false;
      }
    }

    return true;
  }

  /// True if all elements are the same.
  template<typename A, typename Real, size_t dim>
  CODI_INLINE bool operator==(Direction<Real, dim> const& v, A const& s) {
    return s == v;
  }

  /// Output stream operator.
  template<typename Real, size_t dim>
  std::ostream& operator<<(std::ostream& os, Direction<Real, dim> const& v) {
    os << "{";
    for (size_t i = 0; i < dim; ++i) {
      if (i != 0) {
        os << ", ";
      }
      os << v[i];
    }
    os << "}";

    return os;
  }

#ifndef DOXYGEN_DISABLE
  template<typename _Type>
  struct RealTraits::IsTotalZero<_Type, GradientTraits::EnableIfDirection<_Type>> {
    public:

      using Type = CODI_DD(_Type, TEMPLATE(Direction<double, 1>));
      using Real = typename GradientTraits::Real<Type>;

      static CODI_INLINE bool isTotalZero(Type const& v) {
        for (size_t i = 0; i < GradientTraits::dim<Type>(); ++i) {
          if (!codi::RealTraits::isTotalZero(v[i])) {
            return false;
          }
        }
        return true;
      }
  };

  template<typename _Type>
  struct RealTraits::IsTotalFinite<_Type, GradientTraits::EnableIfDirection<_Type>> {
    public:

      using Type = CODI_DD(_Type, TEMPLATE(Direction<double, 1>));

      static CODI_INLINE bool isTotalFinite(Type const& v) {
        for (size_t i = 0; i < GradientTraits::dim<Type>(); ++i) {
          if (!codi::RealTraits::isTotalFinite(v[i])) {
            return false;
          }
        }
        return true;
      }
  };

  namespace GradientTraits {

    template<typename _Gradient>
    struct TraitsImplementation<_Gradient, EnableIfDirection<_Gradient>> {
      public:

        using Gradient = CODI_DD(_Gradient, TEMPLATE(Direction<double, 1>));
        using Real = typename Gradient::Real;

        static size_t constexpr dim = Gradient::dim;

        CODI_INLINE static Real& at(Gradient& gradient, size_t dim) {
          return gradient[dim];
        }

        CODI_INLINE static Real const& at(Gradient const& gradient, size_t dim) {
          return gradient[dim];
        }
    };
  }
#endif
}
