#pragma once

#include "../aux/macros.h"
#include "../config.h"

/** \copydoc codi::Namespace */
namespace codi {

  template<typename _Type, typename = void>
  struct RealTraits {

      using Type = DECLARE_DEFAULT(_Type, double);

      using PassiveReal = Type;

      static int constexpr MaxDerivativeOrder = 0;

      PassiveReal const& getPassiveValue(Type const& v) {
        return v;
      }
  };

  template<typename Type>
  using PassiveRealType = typename RealTraits<Type>::PassiveReal;

  template<typename Type>
  inline int constexpr MaxDerivativeOrderValue = RealTraits<Type>::MaxDerivativeOrder;

  template<typename Type>
  PassiveRealType<Type> const& getPassiveValue(Type const& v) {
    return RealTraits<Type>::getPassiveValue(v);
  }
}
