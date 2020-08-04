#pragma once

#include <cmath>

#include "../../aux/exceptions.hpp"
#include "../../aux/macros.hpp"
#include "../../config.h"
#include "../../traits/realTraits.hpp"
#include "../binaryExpression.hpp"
#include "../constantExpression.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /*******************************************************************************
   * Section: Builtin binary operators
   *
   * Description: TODO
   *
   */

  template<typename _Real>
  struct Add : public BinaryOperation<_Real> {
    public:

      using Real = CODI_DECLARE_DEFAULT(_Real, double);

      template<typename ArgA, typename ArgB>
      static CODI_INLINE Real primal(ArgA const& argA, ArgB const& argB) {
        return argA + argB;
      }

      template<typename ArgA, typename ArgB>
      static CODI_INLINE PassiveRealType<Real> gradientA(ArgA const& argA, ArgB const& argB, Real const& result) {
        CODI_UNUSED(argA, argB, result);
        return 1.0;
      }

      template<typename ArgA, typename ArgB>
      static CODI_INLINE  PassiveRealType<Real> gradientB(ArgA const& argA, ArgB const& argB, Real const& result) {
        CODI_UNUSED(argA, argB, result);
        return 1.0;
      }
  };
  #define OPERATION_LOGIC Add
  #define FUNCTION operator +
  #include "binaryOverloads.tpp"

  template<typename _Real>
  struct Substract : public BinaryOperation<_Real> {
    public:

      using Real = CODI_DECLARE_DEFAULT(_Real, double);

      template<typename ArgA, typename ArgB>
      static CODI_INLINE Real primal(ArgA const& argA, ArgB const& argB) {
        return argA - argB;
      }

      template<typename ArgA, typename ArgB>
      static CODI_INLINE PassiveRealType<Real> gradientA(ArgA const& argA, ArgB const& argB, Real const& result) {
        CODI_UNUSED(argA, argB, result);
        return 1.0;
      }

      template<typename ArgA, typename ArgB>
      static CODI_INLINE  PassiveRealType<Real> gradientB(ArgA const& argA, ArgB const& argB, Real const& result) {
        CODI_UNUSED(argA, argB, result);
        return -1.0;
      }
  };
  #define OPERATION_LOGIC Substract
  #define FUNCTION operator -
  #include "binaryOverloads.tpp"

  template<typename _Real>
  struct Multiply : public BinaryOperation<_Real> {
    public:

      using Real = CODI_DECLARE_DEFAULT(_Real, double);

      template<typename ArgA, typename ArgB>
      static CODI_INLINE Real primal(ArgA const& argA, ArgB const& argB) {
        return argA * argB;
      }

      template<typename ArgA, typename ArgB>
      static CODI_INLINE ArgB const& gradientA(ArgA const& argA, ArgB const& argB, Real const& result) {
        CODI_UNUSED(argA, argB, result);
        return argB;
      }

      template<typename ArgA, typename ArgB>
      static CODI_INLINE ArgA const& gradientB(ArgA const& argA, ArgB const& argB, Real const& result) {
        CODI_UNUSED(argA, argB, result);
        return argA;
      }
  };
  #define OPERATION_LOGIC Multiply
  #define FUNCTION operator *
  #include "binaryOverloads.tpp"

  template<typename _Real>
  struct Divide : public BinaryOperation<_Real> {
    public:

      using Real = CODI_DECLARE_DEFAULT(_Real, double);

      template<typename ArgA, typename ArgB>
      static CODI_INLINE Real primal(ArgA const& argA, ArgB const& argB) {
        return argA / argB;
      }

      template<typename ArgA, typename ArgB>
      static CODI_INLINE Real gradientA(ArgA const& argA, ArgB const& argB, Real const& result) {
        CODI_UNUSED(argA, result);

        checkArguments(argB);
        return 1.0 / argB;
      }

      template<typename ArgA, typename ArgB>
      static CODI_INLINE  Real gradientB(ArgA const& argA, ArgB const& argB, Real const& result) {
        CODI_UNUSED(argA, result);

        checkArguments(argB);
        return -result / argB;
      }

    private:
      template<typename ArgB>
      static CODI_INLINE void checkArguments(ArgB const& argB) {
        if(Config::CheckExpressionArguments) {
          if( isTotalZero(getPassiveValue(argB))) {
            CODI_EXCEPTION("Division called with divisor of zero.");
          }
        }
      }
  };
  #define OPERATION_LOGIC Divide
  #define FUNCTION operator /
  #include "binaryOverloads.tpp"

  /*******************************************************************************
   * Section: Standard math library binary operators
   *
   * Description: TODO
   *
   */

  using std::atan2;
  using std::copysign;
  using std::max;
  using std::min;
  using std::pow;
  using std::remainder;

  template<typename _Real>
  struct Atan2 : public BinaryOperation<_Real> {
    public:

      using Real = CODI_DECLARE_DEFAULT(_Real, double);

      template<typename ArgA, typename ArgB>
      static CODI_INLINE Real primal(ArgA const& argA, ArgB const& argB) {
        return atan2(argA, argB);
      }

      template<typename ArgA, typename ArgB>
      static CODI_INLINE Real gradientA(ArgA const& argA, ArgB const& argB, Real const& result) {
        CODI_UNUSED(result);

        checkArguments(argA, argB);
        Real divisor = argA * argA + argB * argB;
        divisor = 1.0 / divisor;
        return argB * divisor;
      }

      template<typename ArgA, typename ArgB>
      static CODI_INLINE  Real gradientB(ArgA const& argA, ArgB const& argB, Real const& result) {
        CODI_UNUSED(result);

        checkArguments(argA, argB);
        Real divisor = argA * argA + argB * argB;
        divisor = 1.0 / divisor;
        return -argA * divisor;
      }

    private:
      template<typename ArgA, typename ArgB>
      static CODI_INLINE void checkArguments(ArgA& argA, ArgB& argB) {
        if(Config::CheckExpressionArguments) {
          if( 0.0 == getPassiveValue(argA) && 0.0 == getPassiveValue(argB)) {
            CODI_EXCEPTION("atan2 called at point (0,0).");
          }
        }
      }
  };
  #define OPERATION_LOGIC Atan2
  #define FUNCTION atan2
  #include "binaryOverloads.tpp"

  template<typename _Real>
  struct Copysign : public BinaryOperation<_Real> {
    public:

      using Real = CODI_DECLARE_DEFAULT(_Real, double);

      template<typename ArgA, typename ArgB>
      static CODI_INLINE Real primal(ArgA const& argA, ArgB const& argB) {
        return copysign(argA, argB);
      }

      template<typename ArgA, typename ArgB>
      static CODI_INLINE PassiveRealType<Real> gradientA(ArgA const& argA, ArgB const& argB, Real const& result) {

        CODI_UNUSED(result);
        if(getPassiveValue(argA) < 0.0) {
          if (getPassiveValue(argB) < 0.0) {
            return PassiveRealType<Real>(1.0);
          } else {
            return PassiveRealType<Real>(-1.0);
          }
        } else if(getPassiveValue(argA) > 0.0) {
          if (getPassiveValue(argB) < 0.0) {
            return PassiveRealType<Real>(-1.0);
          } else {
            return PassiveRealType<Real>(1.0);
          }
        } else {
          return PassiveRealType<Real>(0.0);
        }
      }

      template<typename ArgA, typename ArgB>
      static CODI_INLINE  PassiveRealType<Real> gradientB(ArgA const& argA, ArgB const& argB, Real const& result) {
        CODI_UNUSED(argA, argB, result);

        return PassiveRealType<Real>(0.0);
      }
  };

  #define OPERATION_LOGIC Copysign
  #define FUNCTION copysign
  #include "binaryOverloads.tpp"

  #define OPERATION_LOGIC Copysign
  #define FUNCTION copysignf
  #include "binaryOverloads.tpp"

  template<typename _Real>
  struct Max : public BinaryOperation<_Real> {
    public:

      using Real = CODI_DECLARE_DEFAULT(_Real, double);

      template<typename ArgA, typename ArgB>
      static CODI_INLINE Real primal(ArgA const& argA, ArgB const& argB) {
        return max(argA, argB);
      }

      template<typename ArgA, typename ArgB>
      static CODI_INLINE PassiveRealType<Real> gradientA(ArgA const& argA, ArgB const& argB, Real const& result) {
        CODI_UNUSED(result);

        if(getPassiveValue(argA) > getPassiveValue(argB)) {
          return 1.0;
        } else {
          return 0.0;
        }
      }

      template<typename ArgA, typename ArgB>
      static CODI_INLINE  PassiveRealType<Real> gradientB(ArgA const& argA, ArgB const& argB, Real const& result) {
        CODI_UNUSED(result);

        if(getPassiveValue(argA) > getPassiveValue(argB)) {
          return 0.0;
        } else {
          return 1.0;
        }
      }
  };

  #define OPERATION_LOGIC Max
  #define FUNCTION max
  #include "binaryOverloads.tpp"

  #define OPERATION_LOGIC Max
  #define FUNCTION fmax
  #include "binaryOverloads.tpp"

  template<typename _Real>
  struct Min : public BinaryOperation<_Real> {
    public:

      using Real = CODI_DECLARE_DEFAULT(_Real, double);

      template<typename ArgA, typename ArgB>
      static CODI_INLINE Real primal(ArgA const& argA, ArgB const& argB) {
        return min(argA, argB);
      }

      template<typename ArgA, typename ArgB>
      static CODI_INLINE PassiveRealType<Real> gradientA(ArgA const& argA, ArgB const& argB, Real const& result) {
        CODI_UNUSED(result);

        if(getPassiveValue(argA) < getPassiveValue(argB)) {
          return 1.0;
        } else {
          return 0.0;
        }
      }

      template<typename ArgA, typename ArgB>
      static CODI_INLINE  PassiveRealType<Real> gradientB(ArgA const& argA, ArgB const& argB, Real const& result) {
        CODI_UNUSED(result);

        if(getPassiveValue(argA) < getPassiveValue(argB)) {
          return 0.0;
        } else {
          return 1.0;
        }
      }
  };
  #define OPERATION_LOGIC Min
  #define FUNCTION min
  #include "binaryOverloads.tpp"

  #define OPERATION_LOGIC Min
  #define FUNCTION fmin
  #include "binaryOverloads.tpp"

  template<typename _Real>
  struct Pow : public BinaryOperation<_Real> {
    public:

      using Real = CODI_DECLARE_DEFAULT(_Real, double);

      template<typename ArgA, typename ArgB>
      static CODI_INLINE Real primal(ArgA const& argA, ArgB const& argB) {
        return pow(argA, argB);
      }

      template<typename ArgA, typename ArgB>
      static CODI_INLINE Real gradientA(ArgA const& argA, ArgB const& argB, Real const& result) {
        CODI_UNUSED(result);

        checkArguments(argA);
        if (getPassiveValue(argA) <= 0.0 && 1 <= RealTraits<ArgB>::MaxDerivativeOrder) {
          // Special case for higher order derivatives. Derivative will be wrong since the argB part is not evaluated.
          return getPassiveValue(argB) * pow(argA, argB - 1.0);
        } else {
          return argB * pow(argA, argB - 1.0);
        }
      }

      template<typename ArgA, typename ArgB>
      static CODI_INLINE  Real gradientB(ArgA const& argA, ArgB const& argB, Real const& result) {
        CODI_UNUSED(argB);

        checkArguments(argA);
        if (getPassiveValue(argA) > 0.0) {
          return log(argA) * result;
        } else {
          return Real();
        }
      }

    private:
      template<typename ArgA>
      static CODI_INLINE void checkArguments(ArgA& argA) {
        if(Config::CheckExpressionArguments) {
          if(getPassiveValue(argA) < 0.0) {
            CODI_EXCEPTION("Negative base for active exponent in pow function. (Value: %0.15e)", getPassiveValue(argA));
          }
        }
      }
  };
  #define OPERATION_LOGIC Pow
  #define FUNCTION pow
  #include "binaryOverloads.tpp"

  /* TODO: Refer to underlying IEEE remainder operation.
   * Denote primal formula.
   */
  template<typename _Real>
  struct Remainder : public BinaryOperation<_Real> {
    public:

      using Real = CODI_DECLARE_DEFAULT(_Real, double);

      template<typename ArgA, typename ArgB>
      static CODI_INLINE Real primal(ArgA const& argA, ArgB const& argB) {
        return remainder(argA, argB);
      }

      template<typename ArgA, typename ArgB>
      static CODI_INLINE Real gradientA(ArgA const& argA, ArgB const& argB, Real const& result) {
        CODI_UNUSED(argA, argB, result);

        return (Real)1.0;
      }

      template<typename ArgA, typename ArgB>
      static CODI_INLINE  Real gradientB(ArgA const& argA, ArgB const& argB, Real const& result) {
        CODI_UNUSED(result);

        checkArguments(argB);

        using std::round;
        return -round(argA/argB);
      }

    private:
      template<typename ArgB>
      static CODI_INLINE void checkArguments(ArgB const& argB) {
        if (Config::CheckExpressionArguments) {
          if (0.0 == getPassiveValue(argB)) {
            CODI_EXCEPTION("Remainder called with divisor of zero.");
          }
        }
      }
  };
  #define OPERATION_LOGIC Remainder
  #define FUNCTION remainder
  #include "binaryOverloads.tpp"
}

namespace std {
  using codi::atan2;
  using codi::copysign;
  using codi::copysignf;
  using codi::fmax;
  using codi::fmin;
  using codi::max;
  using codi::min;
  using codi::pow;
  using codi::remainder;
}
