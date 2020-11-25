#pragma once

#include <cmath>

#include "../aux/macros.hpp"
#include "../config.h"

/** \copydoc codi::Namespace */
namespace codi {

  namespace RealTraits {

    template<typename _Type, typename = void>
    struct TraitsImplementation {
      public:

        using Type = CODI_DECLARE_DEFAULT(_Type, double);

        using PassiveReal = Type;

        static int constexpr MaxDerivativeOrder = 0;

        static CODI_INLINE PassiveReal const& getPassiveValue(Type const& v) {
          return v;
        }
    };

    /// TODO
    template<typename _Type, typename = void>
    struct IsTotalZero {
      public:

        using Type = CODI_DECLARE_DEFAULT(_Type, double);

        /// TODO
        static CODI_INLINE bool isTotalZero(Type const& v) {
          return Type() == v;
        }
    };

    /// TODO
    template<typename _Type, typename = void>
    struct IsTotalFinite {
      public:

        using Type = CODI_DECLARE_DEFAULT(_Type, double);

        /// TODO
        static CODI_INLINE bool isTotalFinite(Type const& v) {
          using std::isfinite;
          return isfinite(v);
        }
    };

    template<typename Type>
    using PassiveReal = typename TraitsImplementation<Type>::PassiveReal;

    template<typename Type>
    CODI_INLINE size_t constexpr MaxDerivativeOrder() {
      return TraitsImplementation<Type>::MaxDerivativeOrder;
    }

    template<typename Type>
    CODI_INLINE PassiveReal<Type> const& getPassiveValue(Type const& v) {
      return TraitsImplementation<Type>::getPassiveValue(v);
    }

    template<typename Type>
    CODI_INLINE bool isTotalZero(Type const& v) {
      return IsTotalZero<Type>::isTotalZero(v);
    }

    template<typename Type>
    CODI_INLINE bool isTotalFinite(Type const& v) {
      return IsTotalFinite<Type>::isTotalFinite(v);
    }

    template<typename Type>
    using IsPassiveReal = std::is_same<Type, PassiveReal<Type>>;

    template<typename Type>
    using enableIfIsNotPassiveReal = typename std::enable_if<!IsPassiveReal<Type>::value>::type;

    template<typename Type>
    using enableIfIsPassiveReal = typename std::enable_if<IsPassiveReal<Type>::value>::type;

  }
}
