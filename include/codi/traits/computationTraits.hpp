#pragma once

#include <complex>
#include <type_traits>

#include "../aux/macros.hpp"
#include "../config.h"

/** \copydoc codi::Namespace */
namespace codi {

  namespace ComputationTraits {

    template<typename _Lhs, typename _Rhs, typename = void>
    struct Update {
      public:
        using Lhs = CODI_DD(_Lhs, double);
        using Rhs = CODI_DD(_Rhs, double);

        using Return = Lhs&;

        static Return update(Lhs& lhs, Rhs const& rhs) {
          return lhs += rhs;
        }
    };

    template<typename Lhs, typename Rhs>
    CODI_INLINE typename Update<Lhs, Rhs>::Return update(Lhs& lhs, Rhs const& rhs) {
      return Update<Lhs, Rhs>::update(lhs, rhs);
    }

    template<typename _Jacobian, typename = void>
    struct Transpose {
      public:
        static_assert(false && std::is_void<_Jacobian>::value, "Instantiation of unspecialized Jacobian transpose.");
        using Jacobian = CODI_DD(_Jacobian, CODI_ANY);

        using Return = CODI_ANY;

        static Return transpose(Jacobian const& jacobian);
    };

    template<typename Jacobian>
    CODI_INLINE typename Transpose<Jacobian>::Return transpose(Jacobian const& jacobian) {
      return Transpose<Jacobian>::transpose(jacobian);
    }
  }

#define CODI_CREATE_UPDATE(LhsType, RhsType, up)       \
  template<>                                           \
  struct ComputationTraits::Update<LhsType, RhsType> { \
    public:                                            \
      using Lhs = LhsType;                             \
      using Rhs = RhsType;                             \
      using Return = LhsType&;                         \
      static Return update(Lhs& lhs, Rhs const& rhs) { \
        return up;                                     \
      }                                                \
  }

#define CODI_CREATE_TRANSPOSE(Type, RType, trans)     \
  template<>                                          \
  struct ComputationTraits::Transpose<Type> {         \
    public:                                           \
      using Jacobian = Type;                          \
      using Return = RType;                           \
      static Return transpose(Type const& jacobian) { \
        return trans;                                 \
      }                                               \
  }

  template<typename T>
  struct ComputationTraits::Transpose<T, typename std::enable_if<std::is_floating_point<T>::value>::type> {
    public:
      using Jacobian = T;
      using Return = T;

      static Return transpose(Jacobian const& jacobian) {
        return jacobian;
      }
  };

  template<typename Inner>
  struct ComputationTraits::Transpose<std::complex<Inner>> {
    public:
      using Jacobian = std::complex<Inner>;
      using Return = std::complex<Inner>;

      static Return transpose(Jacobian const& jacobian) {
        return std::conj(jacobian);
      }
  };

  template<typename Inner>
  struct ComputationTraits::Update<Inner, std::complex<Inner>> {
    public:
      using Lhs = Inner;
      using Rhs = std::complex<Inner>;

      using Return = Inner&;

      static Return update(Lhs& lhs, Rhs const& rhs) {
        return lhs += std::real(rhs);
      }
  };
}
