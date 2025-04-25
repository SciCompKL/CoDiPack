/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2025 Chair for Scientific Computing (SciComp), University of Kaiserslautern-Landau
 * Homepage: http://scicomp.rptu.de
 * Contact:  Prof. Nicolas R. Gauger (codi@scicomp.uni-kl.de)
 *
 * Lead developers: Max Sagebaum, Johannes Blühdorn (SciComp, University of Kaiserslautern-Landau)
 *
 * This file is part of CoDiPack (http://scicomp.rptu.de/software/codi).
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
#include "../misc/self.hpp"
#include "expressionTraits.hpp"
#include "../expressions/complex/complexPredef.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  namespace RealTraits {
    template<typename T_Type, typename>
    struct TraitsImplementation;
  }

  namespace ComputationTraits {

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

/// Create a specialization of ComputationTraits::TransposeImpl
///
/// Has to be used in the codi namespace.
#define CODI_CREATE_TRANSPOSE(Type, RType, trans)     \
  template<>                                          \
  struct ComputationTraits::TransposeImpl<Type> {     \
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

  /// Transpose specialization for expressions of floating point numbers.
  template<typename T>
  struct ComputationTraits::TransposeImpl<
      T, std::enable_if_t<
             ExpressionTraits::isExpression<T> &
             std::is_floating_point<typename RealTraits::TraitsImplementation<T, void>::PassiveReal>::value>> {

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

  /// Transpose specialization for active complex types.
  template<typename Inner>
  struct ComputationTraits::TransposeImpl<ActiveComplex<Inner>> {
    public:
      using Jacobian = ActiveComplex<Inner>;
      using Return = ActiveComplex<Inner>;

      static Return transpose(Jacobian const& jacobian) {
        return conj(jacobian);
      }
  };

  CODI_CREATE_TRANSPOSE(int, int, jacobian);

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

  ///
  /// This is the adjoint of std::complex<double> = double;
  template<typename Inner>
  struct ComputationTraits::UpdateImpl<Inner, ActiveComplex<Inner>> {
    public:
      using Lhs = Inner;
      using Rhs = ActiveComplex<Inner>;

      using Return = Inner&;

      static Return update(Lhs& lhs, Rhs const& rhs) {
        return lhs += real(rhs);
      }
  };
#endif
}
