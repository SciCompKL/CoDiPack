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
    /// Simple specializations can be created with the macro CODI_CREATE_UPDATE.
    ///
    /// @tparam _Lhs Type of the left hand side.
    /// @tparam _Rhs Type of the right hand side.
    template<typename _Lhs, typename _Rhs, typename = void>
    struct UpdateImpl {
      public:
        using Lhs = CODI_DD(_Lhs, double);  ///< See UpdateImpl.
        using Rhs = CODI_DD(_Rhs, double);  ///< See UpdateImpl.

        using Return = Lhs&;  ///< Deduced return type.

        /// Perform lhs += rhs.
        static Return update(Lhs& lhs, Rhs const& rhs) {
          return lhs += rhs;
        }
    };

    /// Perform lhs += rhs. Logic can be specialized by specializing UpdateImpl.
    template<typename Lhs, typename Rhs>
    CODI_INLINE typename UpdateImpl<Lhs, Rhs>::Return update(Lhs& lhs, Rhs const& rhs) {
      return UpdateImpl<Lhs, Rhs>::update(lhs, rhs);
    }

    /// @brief Perform \f$ a^T \f$ or \f$ a^H \f$ if entries are complex. No default implementation available.
    ///
    /// Simple specializations can be created with the macro CODI_CREATE_TRANSPOSE.
    ///
    /// @tparam _Jacobian  Deduced type of the Jacobian.
    template<typename _Jacobian, typename = void>
    struct TransposeImpl {
      public:
        static_assert(false && std::is_void<_Jacobian>::value, "Instantiation of unspecialized Jacobian transpose.");
        using Jacobian = CODI_DD(_Jacobian, CODI_ANY);  ///< See TransposeImpl.

        using Return = CODI_ANY;  ///< Deduced return type.

        /// @brief Perform \f$ a^T \f$ or \f$ a^H \f$ if entries are complex.
        static Return transpose(Jacobian const& jacobian);
    };

    /// Perform \f$ a^T \f$ or \f$ a^H \f$ if entries are complex. Logic can be specialized by specializing
    /// TransposeImpl.
    template<typename Jacobian>
    CODI_INLINE typename TransposeImpl<Jacobian>::Return transpose(Jacobian const& jacobian) {
      return TransposeImpl<Jacobian>::transpose(jacobian);
    }

    /// @brief Perform the adjoint of Outer(Inner). Default implementation for Outer == Inner.
    ///
    /// Simple specialization can be create with the macro CODI_CREATE_ADJOINT_CONVERSION.
    ///
    /// @tparam _Outer  Type into which was converted.
    /// @tparam _Inner  Type from which was converted.
    template<typename _Outer, typename _Inner, typename = void>
    struct AdjointConversionImpl {
      public:
        static_assert(false && std::is_void<_Outer>::value, "Instantiation of unspecialized adjoint conversion.");
        using Outer = CODI_DD(_Outer, CODI_ANY);  ///< See ReduceImpl.
        using Inner = CODI_DD(_Inner, CODI_ANY);  ///< See ReduceImpl.

        using Return = CODI_ANY;  ///< Deduced return type.

        /// @brief Perform the adjoint of Outer(Inner).
        static Return adjointConversion(Outer const& jacobian);
    };

    /// Perform the adjoint of Outer(Inner). Logic can be specialized via AdjointConversionImpl.
    template<typename Inner, typename Outer>
    CODI_INLINE typename AdjointConversionImpl<Outer, Inner>::Return adjointConversion(Outer const& jacobian) {
      return AdjointConversionImpl<Outer, Inner>::adjointConversion(jacobian);
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

/// Create a specialization of ComputationTraits::AdjointConversionImpl
///
/// Has to be used in the codi namespace.
#define CODI_CREATE_ADJOINT_CONVERSION(OuterType, InnerType, RType, conversion) \
  template<>                                                                    \
  struct ComputationTraits::AdjointConversionImpl<OuterType, InnerType> {       \
    public:                                                                     \
      using Outer = OuteType;                                                   \
      using Inner = InnerType;                                                  \
      using Return = RType;                                                     \
      static Return adjointConversion(Outer const& jacobian) {                  \
        return conversion;                                                      \
      }                                                                         \
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

  /// Adjoint conversion specialization for Inner == Outer
  template<typename Type>
  struct ComputationTraits::AdjointConversionImpl<Type, Type> {
    public:
      using Outer = Type;
      using Inner = Type;

      using Return = Type const&;

      static Return adjointConversion(Outer const& jacobian) {
        return jacobian;
      }
  };

  /// Adjoint conversion specialization for Outer == std::complex<Inner>
  template<typename Type>
  struct ComputationTraits::AdjointConversionImpl<std::complex<Type>, Type> {
    public:
      using Outer = std::complex<Type>;
      using Inner = Type;

      using Return = Type;

      static Return adjointConversion(Outer const& jacobian) {
        return std::real(jacobian);
      }
  };
#endif
}
