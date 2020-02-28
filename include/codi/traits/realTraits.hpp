#pragma once

#include "../aux/macros.h"
#include "../config.h"

/** \copydoc codi::Namespace */
namespace codi {

  template<typename _Type, typename = void>
  struct RealTraits {
    public:

      using Type = DECLARE_DEFAULT(_Type, double);

      using PassiveReal = Type;

      static int constexpr MaxDerivativeOrder = 0;

      static CODI_INLINE PassiveReal const& getPassiveValue(Type const& v) {
        return v;
      }

      static CODI_INLINE bool isTotalZero(Type const& v) {
        return Type() == v;
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
    return RealTraits<Type>::isTotalZero(v);
  }
}
