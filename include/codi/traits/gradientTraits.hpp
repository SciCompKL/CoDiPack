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

#include <array>
#include <type_traits>

#include "../misc/macros.hpp"
#include "atomicTraits.hpp"
#include "misc/enableIfHelpers.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  template<typename Real, size_t dim>
  struct Direction;

  /// Traits for everything that can be an used as a gradient (adjoint, tangent) usually the second template argument
  /// of CoDi Gen types codi::RealReverseGen. Possible types are double, codi::RealReverse, codi::Direction etc..
  namespace GradientTraits {

    /*******************************************************************************/
    /// @name General gradient traits
    /// @{

    /**
     * @brief Common traits for all types used as gradients
     *
     * @tparam T_Gradient  The type of the gradient.
     */
    template<typename T_Gradient, typename = void>
    struct TraitsImplementation {
      public:

        using Gradient = CODI_DD(T_Gradient, double);  ///< See TraitsImplementation
        using Real = Gradient;                         ///< The base value used in the gradient entries.

        static size_t constexpr dim = 1;  ///< Number of dimensions this gradient value has.

        /// Get the entry at the given index.
        CODI_INLINE static Real& at(Gradient& gradient, size_t dim) {
          CODI_UNUSED(dim);
          return gradient;
        }

        /// Get the entry at the given index.
        CODI_INLINE static Real const& at(Gradient const& gradient, size_t dim) {
          CODI_UNUSED(dim);
          return gradient;
        }

        /// Converts the (possibly multi-component) gradient to an array of Reals.
        CODI_INLINE static std::array<AtomicTraits::RemoveAtomic<Real>, dim> toArray(Gradient const& gradient) {
          return std::array<AtomicTraits::RemoveAtomic<Real>, dim>{at(gradient, 0)};
        }
    };

    /// \copydoc codi::GradientTraits::TraitsImplementation::Real
    template<typename Gradient>
    using Real = typename TraitsImplementation<Gradient>::Real;

    /// \copydoc codi::GradientTraits::TraitsImplementation::dim
    template<typename Gradient>
    CODI_INLINE size_t constexpr dim() {
      return TraitsImplementation<Gradient>::dim;
    }

    /// \copydoc codi::GradientTraits::TraitsImplementation::at()
    template<typename Gradient>
    CODI_INLINE typename TraitsImplementation<Gradient>::Real& at(Gradient& gradient, size_t dim) {
      return TraitsImplementation<Gradient>::at(gradient, dim);
    }

    /// \copydoc codi::GradientTraits::TraitsImplementation::at()
    template<typename Gradient>
    CODI_INLINE typename TraitsImplementation<Gradient>::Real const& at(Gradient const& gradient, size_t dim) {
      return TraitsImplementation<Gradient>::at(gradient, dim);
    }

    /// \copydoc codi::GradientTraits::TraitsImplementation::toArray()
    template<typename Gradient>
    CODI_INLINE std::array<AtomicTraits::RemoveAtomic<typename TraitsImplementation<Gradient>::Real>,
                           TraitsImplementation<Gradient>::dim>
    toArray(Gradient const& gradient) {
      return TraitsImplementation<Gradient>::toArray(gradient);
    }

    /// @}
    /*******************************************************************************/
    /// @name Detection of specific gradient types
    /// @{

    /// If the expression inherits from Direction. Is either std::false_type or std::true_type
    template<typename Gradient, typename = void>
    struct IsDirection : std::false_type {};

#ifndef DOXYGEN_DISABLE
    template<typename Gradient>
    struct IsDirection<Gradient,
                       typename enable_if_same<Gradient, Direction<typename Gradient::Real, Gradient::dim> >::type>
        : std::true_type {};
#endif

#if CODI_IS_CPP14
    /// Value entry of IsDirection
    template<typename Gradient>
    bool constexpr isDirection = IsDirection<Gradient>::value;
#endif

    /// Enable if wrapper for EnableIfDirection
    template<typename Gradient>
    using EnableIfDirection = typename std::enable_if<IsDirection<Gradient>::value>::type;

    /// @}
  }
}
