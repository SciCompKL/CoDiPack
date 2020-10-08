#pragma once

#include <type_traits>

#include "../aux/macros.hpp"
#include "aux/enableIfHelpers.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  template<typename Real, size_t dim>
  struct Direction;

  namespace GradientTraits {

    template<typename _Gradient, typename = void>
    struct TraitsImplementation {
      public:

        using Gradient = CODI_DECLARE_DEFAULT(_Gradient, double);
        using Real = Gradient;

        static size_t constexpr dim = 1;

        CODI_INLINE static Real& at(Gradient& gradient, size_t dim) {
          CODI_UNUSED(dim);
          return gradient;
        }

        CODI_INLINE static Real const& at(Gradient const& gradient, size_t dim) {
          CODI_UNUSED(dim);
          return gradient;
        }
    };

    template<typename Gradient>
    using Real = typename TraitsImplementation<Gradient>::Real;

    template<typename Gradient>
    CODI_INLINE size_t constexpr dim() {
      return TraitsImplementation<Gradient>::dim;
    }

    template<typename Gradient>
    CODI_INLINE typename TraitsImplementation<Gradient>::Real& at(Gradient& gradient, size_t dim) {
      return TraitsImplementation<Gradient>::at(gradient, dim);
    }

    template<typename Gradient>
    CODI_INLINE typename TraitsImplementation<Gradient>::Real const& at(Gradient const& gradient, size_t dim) {
      return TraitsImplementation<Gradient>::at(gradient, dim);
    }

    template<typename Gradient, typename = void>
    struct IsDirection : std::false_type {};

    template<typename Gradient>
    struct IsDirection<
      Gradient,
      typename enable_if_same<
        Gradient,
        Direction<typename Gradient::Real, Gradient::dim>
      >::type
    > : std::true_type {};

    template<typename Gradient>
    using isDirection = IsDirection<Gradient>;

    template<typename Gradient>
    using enableIfDirection = typename std::enable_if<isDirection<Gradient>::value>::type;
  }
}
