#pragma once

#include <complex>
#include <type_traits>

#include "../aux/macros.hpp"
#include "../config.h"

/** \copydoc codi::Namespace */
namespace codi {

  namespace ComputationTraits {

    /// @brief Perform the operation lhs += rhs. Default logic uses operator +=.
    ///
    /// Simple specialization can be create with the macro CODI_CREATE_UPDATE.
    ///
    /// @tparam _Lhs Type of the left hand side.
    /// @tparam _Rhs Type of the right hand side.
    template<typename _Lhs, typename _Rhs, typename = void>
    struct UpdateImpl {
      public:
        using Lhs = CODI_DD(_Lhs, double);  ///< See UpdateImpl.
        using Rhs = CODI_DD(_Rhs, double);  ///< See UpdateImpl.

        using Return = Lhs&;  ///< Deduced return type.

        /// Perform lhs += rhs
        static Return update(Lhs& lhs, Rhs const& rhs) {
          return lhs += rhs;
        }
    };

    /// Perform lhs += rhs. Logic can be specialized by specializing UpdateImpl.
    template<typename Lhs, typename Rhs>
    CODI_INLINE typename UpdateImpl<Lhs, Rhs>::Return update(Lhs& lhs, Rhs const& rhs) {
      return UpdateImpl<Lhs, Rhs>::update(lhs, rhs);
    }

    /// @brief Perform \f$ a^T \f$. No default implementation available.
    ///
    /// Simple specialization can be create with the macro CODI_CREATE_TRANSPOSE.
    ///
    /// @tparam _Jacobian  Deduced type of the Jacobian.
    template<typename _Jacobian, typename = void>
    struct TransposeImpl {
      public:
        static_assert(false && std::is_void<_Jacobian>::value, "Instantiation of unspecialized Jacobian transpose.");
        using Jacobian = CODI_DD(_Jacobian, CODI_ANY);  ///< See TransposeImpl.

        using Return = CODI_ANY;  ///< Deduced return type.

        /// @brief Perform \f$ a^T \f$.
        static Return transpose(Jacobian const& jacobian);
    };

    /// Perform \f$ a^T \f$. Logic can be specialized by specializing TransposeImpl.
    template<typename Jacobian>
    CODI_INLINE typename TransposeImpl<Jacobian>::Return transpose(Jacobian const& jacobian) {
      return TransposeImpl<Jacobian>::transpose(jacobian);
    }
  }

/// Create a specialization of ComputationTraits::UpdateImpl
///
/// Has to be used in the codi namespace.
#define CODI_CREATE_UPDATE(LhsType, RhsType, up)           \
  template<>                                               \
  struct ComputationTraits::UpdateImpl<LhsType, RhsType> { \
    public:                                                \
      using Lhs = LhsType;                                 \
      using Rhs = RhsType;                                 \
      using Return = LhsType&;                             \
      static Return update(Lhs& lhs, Rhs const& rhs) {     \
        return up;                                         \
      }                                                    \
  }

/// Create a specialization of ComputationTraits::TransposeImpl
///
/// Has to be used in the codi namespace.
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

#ifndef DOXYGEN_DISABLE

  /// Transpose specialization for floating point numbers.
  template<typename T>
  struct ComputationTraits::TransposeImpl<T, typename std::enable_if<std::is_floating_point<T>::value>::type> {
    public:
      using Jacobian = T;
      using Return = T;

      static Return transpose(Jacobian const& jacobian) {
        return jacobian;
      }
  };

  /// Transpose specialization for complex types.
  template<typename Inner>
  struct ComputationTraits::TransposeImpl<std::complex<Inner>> {
    public:
      using Jacobian = std::complex<Inner>;
      using Return = std::complex<Inner>;

      static Return transpose(Jacobian const& jacobian) {
        return std::conj(jacobian);
      }
  };

  /// Update specialization for double += std::complex<double>
  ///
  /// This is the adjoint of std::complex<double> = double;
  template<typename Inner>
  struct ComputationTraits::UpdateImpl<Inner, std::complex<Inner>> {
    public:
      using Lhs = Inner;
      using Rhs = std::complex<Inner>;

      using Return = Inner&;

      static Return update(Lhs& lhs, Rhs const& rhs) {
        return lhs += std::real(rhs);
      }
  };
#endif
}
