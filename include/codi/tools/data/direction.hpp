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

#include <initializer_list>

#include "../../config.h"
#include "../../misc/macros.hpp"
#include "../../traits/atomicTraits.hpp"
#include "../../traits/gradientTraits.hpp"
#include "../../traits/realTraits.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Fixed size vector mode implementation.
   *
   * Can be used as the gradient template argument in active CoDiPack types.
   *
   * @tparam T_Real  Type of the vector entries.
   * @tparam T_dim  Dimension of the vector mode.
   */
  template<typename T_Real, size_t T_dim>
  struct Direction {
    public:

      using Real = CODI_DD(T_Real, double);  ///< See Direction.

      static size_t constexpr dim = T_dim;  ///< See Direction.

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
        Real const* array = l.begin();  // This is possible because the standard requires an array storage.
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

      /// Assignment operator.
      CODI_INLINE Direction& operator=(Direction const& v) {
        for (size_t i = 0; i < dim; ++i) {
          this->vector[i] = v.vector[i];
        }

        return *this;
      }

      /// Update operator.
      CODI_INLINE Direction& operator+=(Direction const& v) {
        for (size_t i = 0; i < dim; ++i) {
          this->vector[i] += v.vector[i];
        }

        return *this;
      }

      /// Update operator.
      CODI_INLINE Direction& operator-=(Direction const& v) {
        for (size_t i = 0; i < dim; ++i) {
          this->vector[i] -= v.vector[i];
        }

        return *this;
      }
  };

  template<typename Real, size_t dim>
  size_t constexpr Direction<Real, dim>::dim;

  /// Multiplication with a scalar.
  template<typename Real, size_t dim>
  CODI_INLINE Direction<Real, dim> operator*(Real const& s, Direction<Real, dim> const& v) {
    Direction<Real, dim> r;
    for (size_t i = 0; i < dim; ++i) {
      r[i] = s * v[i];
    }

    return r;
  }

  /// Multiplication with a passive scalar.
  template<typename Real, size_t dim, typename = RealTraits::EnableIfNotPassiveReal<Real>>
  CODI_INLINE Direction<Real, dim> operator*(RealTraits::PassiveReal<Real> const& s, Direction<Real, dim> const& v) {
    Direction<Real, dim> r;
    for (size_t i = 0; i < dim; ++i) {
      r[i] = s * v[i];
    }

    return r;
  }

  /// Multiplication of a non-atomic scalar with a direction that has atomic reals.
  template<typename Real, size_t dim, typename = AtomicTraits::EnableIfAtomic<Real>>
  CODI_INLINE Direction<Real, dim> operator*(AtomicTraits::RemoveAtomic<Real> const& s, Direction<Real, dim> const& v) {
    Direction<Real, dim> r;
    for (size_t i = 0; i < dim; ++i) {
      r[i] = s * v[i];
    }

    return r;
  }

  /// Multiplication with a scalar.
  template<typename Real, size_t dim>
  CODI_INLINE Direction<Real, dim> operator*(Direction<Real, dim> const& v, Real const& s) {
    return s * v;
  }

  /// Multiplication with passive a scalar.
  template<typename Real, size_t dim, typename = RealTraits::EnableIfNotPassiveReal<Real>>
  CODI_INLINE Direction<Real, dim> operator*(Direction<Real, dim> const& v, RealTraits::PassiveReal<Real> const& s) {
    return s * v;
  }

  /// Multiplication of a non-atomic scalar with a direction that has atomic reals.
  template<typename Real, size_t dim, typename = AtomicTraits::EnableIfAtomic<Real>>
  CODI_INLINE Direction<Real, dim> operator*(Direction<Real, dim> const& v, AtomicTraits::RemoveAtomic<Real> const& s) {
    return s * v;
  }

  /// Division by a scalar.
  template<typename Real, size_t dim>
  CODI_INLINE Direction<Real, dim> operator/(Direction<Real, dim> const& v, Real const& s) {
    Direction<Real, dim> r;
    for (size_t i = 0; i < dim; ++i) {
      r[i] = v[i] / s;
    }

    return r;
  }

  /// Division by a passive scalar.
  template<typename Real, size_t dim, typename = RealTraits::EnableIfNotPassiveReal<Real>>
  CODI_INLINE Direction<Real, dim> operator/(Direction<Real, dim> const& v, RealTraits::PassiveReal<Real> const& s) {
    Direction<Real, dim> r;
    for (size_t i = 0; i < dim; ++i) {
      r[i] = v[i] / s;
    }

    return r;
  }

  /// Division of a direction with atomic reals by a non-atomic scalar.
  template<typename Real, size_t dim, typename = AtomicTraits::EnableIfAtomic<Real>>
  CODI_INLINE Direction<Real, dim> operator/(Direction<Real, dim> const& v, AtomicTraits::RemoveAtomic<Real> const& s) {
    Direction<Real, dim> r;
    for (size_t i = 0; i < dim; ++i) {
      r[i] = v[i] / s;
    }

    return r;
  }

  /// Summation of two vectors.
  template<typename Real, size_t dim>
  CODI_INLINE Direction<Real, dim> operator+(Direction<Real, dim> const& v1, Direction<Real, dim> const& v2) {
    Direction<Real, dim> r;
    for (size_t i = 0; i < dim; ++i) {
      r[i] = v1[i] + v2[i];
    }

    return r;
  }

  /// Subtraction of two vectors.
  template<typename Real, size_t dim>
  CODI_INLINE Direction<Real, dim> operator-(Direction<Real, dim> const& v1, Direction<Real, dim> const& v2) {
    Direction<Real, dim> r;
    for (size_t i = 0; i < dim; ++i) {
      r[i] = v1[i] - v2[i];
    }

    return r;
  }

  /// Negation.
  template<typename Real, size_t dim>
  CODI_INLINE Direction<Real, dim> operator-(Direction<Real, dim> const& v) {
    Direction<Real, dim> r;
    for (size_t i = 0; i < dim; ++i) {
      r[i] = -v[i];
    }

    return r;
  }

  /// Component-wise test for inequality. True if at least one component differs.
  template<typename Real, size_t dim>
  CODI_INLINE bool operator!=(Direction<Real, dim> const& v1, Direction<Real, dim> const& v2) {
    return !(v1 == v2);
  }

  /// Component-wise test for inequality with scalar. True if at least one component differs.
  template<typename A, typename Real, size_t dim>
  CODI_INLINE bool operator!=(A const& s, Direction<Real, dim> const& v) {
    return !(s == v);
  }

  /// Component-wise test for inequality with scalar. True if at least one component differs.
  template<typename A, typename Real, size_t dim>
  CODI_INLINE bool operator!=(Direction<Real, dim> const& v, A const& s) {
    return s != v;
  }

  /// Component-wise test for equality. True if all components match.
  template<typename Real, size_t dim>
  CODI_INLINE bool operator==(Direction<Real, dim> const& v1, Direction<Real, dim> const& v2) {
    for (size_t i = 0; i < dim; ++i) {
      if (v1[i] != v2[i]) {
        return false;
      }
    }

    return true;
  }

  /// Component-wise test for equality with scalar. True if all components match.
  template<typename A, typename Real, size_t dim>
  CODI_INLINE bool operator==(A const& s, Direction<Real, dim> const& v) {
    for (size_t i = 0; i < dim; ++i) {
      if (s != v[i]) {
        return false;
      }
    }

    return true;
  }

  /// Component-wise test for equality with scalar. True if all components match.
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
  template<typename T_Type>
  struct RealTraits::IsTotalZero<T_Type, GradientTraits::EnableIfDirection<T_Type>> {
    public:

      using Type = CODI_DD(T_Type, CODI_T(Direction<double, 1>));
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

  template<typename T_Type>
  struct RealTraits::IsTotalFinite<T_Type, GradientTraits::EnableIfDirection<T_Type>> {
    public:

      using Type = CODI_DD(T_Type, CODI_T(Direction<double, 1>));

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

    template<typename T_Gradient>
    struct TraitsImplementation<T_Gradient, EnableIfDirection<T_Gradient>> {
      public:

        using Gradient = CODI_DD(T_Gradient, CODI_T(Direction<double, 1>));
        using Real = typename Gradient::Real;

        static size_t constexpr dim = Gradient::dim;

        CODI_INLINE static Real& at(Gradient& gradient, size_t dim) {
          return gradient[dim];
        }

        CODI_INLINE static Real const& at(Gradient const& gradient, size_t dim) {
          return gradient[dim];
        }

        CODI_INLINE static std::array<AtomicTraits::RemoveAtomic<Real>, dim> toArray(Gradient const& gradient) {
          std::array<AtomicTraits::RemoveAtomic<Real>, dim> result;
          for (size_t i = 0; i < dim; ++i) {
            result[i] = at(gradient, i);
          }
          return result;
        }
    };
  }
#endif
}
