#pragma once

#include "../aux/macros.h"
#include "../config.h"

/** \copydoc codi::Namespace */
namespace codi {

  template<typename _Type, typename = void>
  struct RealTraits {

      using Type = DECLARE_DEFAULT(_Type, double);

      using PassiveReal = Type;

      PassiveReal const& getPassiveValue(Type const& v) {
        return v;
      }
  };

  template<typename Type>
  using PassiveRealType = typename RealTraits<Type>::PassiveReal;

  template<typename Type>
  PassiveRealType<Type> const& getPassiveValue(Type const& v) {
    return RealTraits<Type>::getOriginalValue(v);
  }
}
