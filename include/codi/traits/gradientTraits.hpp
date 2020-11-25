#pragma once

#include <type_traits>

#include "../aux/macros.hpp"
#include "aux/enableIfHelpers.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  template<typename Real, size_t dim>
  struct Direction;

  namespace GradientTraits {

    /*******************************************************************************/
    /// @name General gradient traits
    /// @{

    /**
     * @brief Common traits for all types used as gradients
     *
     * @tparam _Gradient  The type of the gradient.
     */
    template<typename _Gradient, typename = void>
    struct TraitsImplementation {
      public:

        using Gradient = CODI_DECLARE_DEFAULT(_Gradient, double); ///< See TraitsImplementation
        using Real = Gradient;  ///< The base value used in the gradient entries.

        static size_t constexpr dim = 1; ///< Number of dimensions this gradient value has.

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

    /// @}
    /*******************************************************************************/
    /// @name Detection of specific gradient types
    /// @{

    /// If the expression inherits from Direction. Is either std::false_type or std::true_type
    template<typename Gradient, typename = void>
    struct IsDirection : std::false_type {};

#ifndef DOXYGEN_DISABLE
    template<typename Gradient>
    struct IsDirection<
      Gradient,
      typename enable_if_same<
        Gradient,
        Direction<typename Gradient::Real, Gradient::dim>
      >::type
    > : std::true_type {};
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
