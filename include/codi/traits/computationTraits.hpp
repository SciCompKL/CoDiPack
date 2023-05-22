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

#include <complex>
#include <type_traits>

#include "../config.h"
#include "../misc/macros.hpp"
#include "expressionTraits.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  namespace ComputationTraits {

    /**
     * @brief Perform the conversion of Outer(Inner) in the adjoint context. Default implementation for Outer == Inner.
     *
     * E.g. for double(std::complex<double>) the conversion returns the real part of the complex number. This is the
     * adjoint of std::complex<double>(double).
     *
     * Simple specialization can be create with the macro CODI_CREATE_ADJOINT_CONVERSION.
     *
     * @tparam T_Outer  Type into which was converted.
     * @tparam T_Inner  Type from which was converted.
     */
    template<typename T_Outer, typename T_Inner, typename = void>
    struct AdjointConversionImpl {
      public:
        CODI_STATIC_ASSERT(false && std::is_void<T_Outer>::value, "Instantiation of unspecialized adjoint conversion.");
        using Outer = CODI_DD(T_Outer, CODI_ANY);  ///< See ReduceImpl.
        using Inner = CODI_DD(T_Inner, CODI_ANY);  ///< See ReduceImpl.

        using Return = CODI_ANY;  ///< Deduced return type.

        /// @brief Perform the adjoint of Outer(Inner).
        static Return adjointConversion(Inner const& jacobian);
    };

    /// Perform the conversion of Outer(Inner) in the adjoint context. Logic can be specialized via
    /// AdjointConversionImpl.
    ///
    /// E.g. for double(std::complex<double>) the conversion returns the real part of the complex number. This is the
    /// adjoint of std::complex<double>(double).
    template<typename Outer, typename Inner>
    CODI_INLINE typename AdjointConversionImpl<Outer, Inner>::Return adjointConversion(Inner const& jacobian) {
      return AdjointConversionImpl<Outer, Inner>::adjointConversion(jacobian);
    }

    /**
     * @brief Perform \f$ a^T \f$ or \f$ a^H \f$ if entries are complex. No default implementation available.
     *
     * Simple specializations can be created with the macro CODI_CREATE_TRANSPOSE.
     *
     * @tparam T_Jacobian  Deduced type of the Jacobian.
     */
    template<typename T_Jacobian, typename = void>
    struct TransposeImpl {
      public:
        CODI_STATIC_ASSERT(false && std::is_void<T_Jacobian>::value,
                           "Instantiation of unspecialized Jacobian transpose.");
        using Jacobian = CODI_DD(T_Jacobian, CODI_ANY);  ///< See TransposeImpl.

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

    /**
     * @brief Perform the operation lhs += rhs. Default logic uses operator +=.
     *
     * Simple specializations can be created with the macro CODI_CREATE_UPDATE.
     *
     * @tparam T_Lhs Type of the left hand side.
     * @tparam T_Rhs Type of the right hand side.
     */
    template<typename T_Lhs, typename T_Rhs, typename = void>
    struct UpdateImpl {
      public:
        using Lhs = CODI_DD(T_Lhs, double);  ///< See UpdateImpl.
        using Rhs = CODI_DD(T_Rhs, double);  ///< See UpdateImpl.

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
      static Return adjointConversion(Inner const& jacobian) {                  \
        return conversion;                                                      \
      }                                                                         \
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

#ifndef DOXYGEN_DISABLE

  /// Adjoint conversion specialization for Inner == Outer
  template<typename Type>
  struct ComputationTraits::AdjointConversionImpl<Type, Type> {
    public:
      using Outer = Type;
      using Inner = Type;

      using Return = Type;

      static Return adjointConversion(Inner const& jacobian) {
        return jacobian;
      }
  };

  /// Adjoint conversion specialization for Inner == codi::Expression
  template<typename T_Outer, typename T_Inner>
  struct ComputationTraits::AdjointConversionImpl<
      T_Outer, T_Inner,
      typename std::enable_if<!std::is_same<T_Outer, T_Inner>::value &
                              ExpressionTraits::IsExpression<T_Inner>::value>::type> {
    public:

      using Outer = T_Outer;
      using Inner = T_Inner;

      using InnerActive = typename Inner::ActiveResult;
      using InnerActiveConversion = ComputationTraits::AdjointConversionImpl<Outer, InnerActive>;

      using Return = typename InnerActiveConversion::Return;

      static Return adjointConversion(Inner const& jacobian) {
        return InnerActiveConversion::adjointConversion(jacobian);
      }
  };

  /// Adjoint conversion specialization for Outer == std::complex<Inner>
  template<typename Type>
  struct ComputationTraits::AdjointConversionImpl<Type, std::complex<Type>> {
    public:
      using Outer = Type;
      using Inner = std::complex<Type>;

      using Return = Type;

      static Return adjointConversion(Inner const& jacobian) {
        return std::real(jacobian);
      }
  };

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
