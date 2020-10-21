#pragma once

#include <cmath>

#include "../aux/macros.hpp"
#include "../config.h"

/** \copydoc codi::Namespace */
namespace codi {

  template<typename _Type, typename = void>
  struct RealTraits {
    public:

      using Type = CODI_DECLARE_DEFAULT(_Type, double);

      using PassiveReal = Type;

      static int constexpr MaxDerivativeOrder = 0;

      static CODI_INLINE PassiveReal const& getPassiveValue(Type const& v) {
        return v;
      }
  };

  template<typename _Type, typename = void>
  struct IsTotalZero {
    public:

      using Type = CODI_DECLARE_DEFAULT(_Type, double);

      static CODI_INLINE bool isTotalZero(Type const& v) {
        return Type() == v;
      }
  };

  template<typename _Type, typename = void>
  struct IsTotalFinite {
    public:

      using Type = CODI_DECLARE_DEFAULT(_Type, double);

      static CODI_INLINE bool isTotalFinite(Type const& v) {
        using std::isfinite;
        return isfinite(v);
      }
  };

  template<typename Type>
  using PassiveRealType = typename RealTraits<Type>::PassiveReal;

  template<typename Type>
  CODI_INLINE PassiveRealType<Type> const& getPassiveValue(Type const& v) {
    return RealTraits<Type>::getPassiveValue(v);
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
  using RealIsPassiveReal = std::is_same<Type, PassiveRealType<Type>>;

  template<typename Type>
  using enableIfRealIsNotPassiveReal = typename std::enable_if<!RealIsPassiveReal<Type>::value>::type;
}
