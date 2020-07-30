#pragma once

#include <cmath>

#include "../../aux/exceptions.hpp"
#include "../../aux/macros.hpp"
#include "../../config.h"
#include "../../traits/realTraits.hpp"
#include "../unaryExpression.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /*******************************************************************************
   * Section: Builtin unary operators
   *
   * Description: TODO
   *
   */

  template<typename _Real>
  struct UnaryMinus : public UnaryOperation<_Real> {
    public:

      using Real = CODI_DECLARE_DEFAULT(_Real, double);

      /** \copydoc UnaryOperation::primal */
      template<typename Arg>
      static CODI_INLINE Real primal(Arg const& arg) {
        return -arg;
      }

      /** \copydoc UnaryOperation::gradient */
      template<typename Arg>
      static CODI_INLINE PassiveRealType<Real> gradient(Arg const& arg, Real const& result) {
        CODI_UNUSED(arg);
        CODI_UNUSED(result);
        return -1.0;
      }
  };
  #define OPERATION_LOGIC UnaryMinus
  #define FUNCTION operator -
  #include "unaryOverloads.tpp"

  template<typename Real, typename Arg>
  CODI_INLINE ExpressionInterface<Real, Arg> const& operator+(ExpressionInterface<Real, Arg> const& arg) {
    return arg;
  }

  /*******************************************************************************
   * Section: Standard math library unary operators
   *
   * Description: TODO
   *
   */

  using std::abs;
  using std::fabs;
  using std::acos;
  using std::asin;
  using std::atan;
  using std::atanh;
  using std::cbrt;
  using std::ceil;
  using std::cos;
  using std::cosh;
  using std::erf;
  using std::erfc;
  using std::exp;
  using std::floor;
  using std::isfinite;
  using std::isinf;
  using std::isnan;
  using std::isnormal;
  using std::log;
  using std::log10;
  using std::sin;
  using std::sinh;
  using std::sqrt;
  using std::tan;
  using std::tanh;
  using std::tgamma;

  template<typename _Real>
  struct Abs : public UnaryOperation<_Real> {
    public:

      using Real = CODI_DECLARE_DEFAULT(_Real, double);

      /** \copydoc UnaryOperation::primal */
      template<typename Arg>
      static CODI_INLINE Real primal(Arg const& arg) {
        return abs(arg);
      }

      /** \copydoc UnaryOperation::gradient */
      template<typename Arg>
      static CODI_INLINE Real gradient(Arg const& arg, Real const& result) {
        CODI_UNUSED(result);
        if(arg < 0.0) {
          return (Real)-1.0;
        } else if(arg > 0.0) {
          return (Real)1.0;
        } else {
          return (Real)0.0;
        }
      }
  };
  #define OPERATION_LOGIC Abs
  #define FUNCTION abs
  #include "unaryOverloads.tpp"

  #define OPERATION_LOGIC Abs
  #define FUNCTION fabs
  #include "unaryOverloads.tpp"

  template<typename _Real>
  struct Acos : public UnaryOperation<_Real> {
    public:

      using Real = CODI_DECLARE_DEFAULT(_Real, double);

      /** \copydoc UnaryOperation::primal */
      template<typename Arg>
      static CODI_INLINE Real primal(Arg const& arg) {
        return acos(arg);
      }

      /** \copydoc UnaryOperation::gradient */
      template<typename Arg>
      static CODI_INLINE Real gradient(Arg const& arg, Real const& result) {
        CODI_UNUSED(result);
        if(Config::CheckExpressionArguments) {
          if(getPassiveValue(arg) <= -1.0 || 1.0 <= getPassiveValue(arg)) {
            CODI_EXCEPTION("acos outside of (-1, 1).(Value: %0.15e)", getPassiveValue(arg));
          }
        }
        return -1.0 / sqrt(1.0 - arg * arg);
      }
  };
  #define OPERATION_LOGIC Acos
  #define FUNCTION acos
  #include "unaryOverloads.tpp"

  template<typename _Real>
  struct Asin : public UnaryOperation<_Real> {
    public:

      using Real = CODI_DECLARE_DEFAULT(_Real, double);

      /** \copydoc UnaryOperation::primal */
      template<typename Arg>
      static CODI_INLINE Real primal(Arg const& arg) {
        return asin(arg);
      }

      /** \copydoc UnaryOperation::gradient */
      template<typename Arg>
      static CODI_INLINE Real gradient(Arg const& arg, Real const& result) {
        CODI_UNUSED(result);
        if(Config::CheckExpressionArguments) {
          if(getPassiveValue(arg) <= -1.0 || 1.0 <= getPassiveValue(arg)) {
            CODI_EXCEPTION("asin outside of (-1, 1).(Value: %0.15e)", getPassiveValue(arg));
          }
        }
        return 1.0 / sqrt(1.0 - arg * arg);
      }
  };
  #define OPERATION_LOGIC Asin
  #define FUNCTION asin
  #include "unaryOverloads.tpp"

  template<typename _Real>
  struct Atan : public UnaryOperation<_Real> {
    public:

      using Real = CODI_DECLARE_DEFAULT(_Real, double);

      /** \copydoc UnaryOperation::primal */
      template<typename Arg>
      static CODI_INLINE Real primal(Arg const& arg) {
        return atan(arg);
      }

      /** \copydoc UnaryOperation::gradient */
      template<typename Arg>
      static CODI_INLINE Real gradient(Arg const& arg, Real const& result) {
        CODI_UNUSED(result);
        return 1.0 / (1 + arg * arg);
      }
  };
  #define OPERATION_LOGIC Atan
  #define FUNCTION atan
  #include "unaryOverloads.tpp"

  template<typename _Real>
  struct Atanh : public UnaryOperation<_Real> {
    public:

      using Real = CODI_DECLARE_DEFAULT(_Real, double);

      /** \copydoc UnaryOperation::primal */
      template<typename Arg>
      static CODI_INLINE Real primal(Arg const& arg) {
        return atanh(arg);
      }

      /** \copydoc UnaryOperation::gradient */
      template<typename Arg>
      static CODI_INLINE Real gradient(Arg const& arg, Real const& result) {
        CODI_UNUSED(result);
        if(Config::CheckExpressionArguments) {
          if(getPassiveValue(arg) <= -1.0 || 1.0 <= getPassiveValue(arg)) {
            CODI_EXCEPTION("atanh outside of (-1, 1).(Value: %0.15e)", getPassiveValue(arg));
          }
        }
        return 1.0 / (1 - arg * arg);
      }
  };
  #define OPERATION_LOGIC Atanh
  #define FUNCTION atanh
  #include "unaryOverloads.tpp"

  template<typename _Real>
  struct Cbrt : public UnaryOperation<_Real> {
    public:

      using Real = CODI_DECLARE_DEFAULT(_Real, double);

      /** \copydoc UnaryOperation::primal */
      template<typename Arg>
      static CODI_INLINE Real primal(Arg const& arg) {
        return cbrt(arg);
      }

      /** \copydoc UnaryOperation::gradient */
      template<typename Arg>
      static CODI_INLINE Real gradient(Arg const& arg, Real const& result) {
        if(Config::CheckExpressionArguments) {
          if(0.0 == getPassiveValue(arg)) {
            CODI_EXCEPTION("Cbrt of zero value.(Value: %0.15e)", getPassiveValue(arg));
          }
        }
        if(result != 0.0) {
          return 1.0 / (3.0 * result * result);
        } else {
          return (Real)0.0;
        }
      }
  };
  #define OPERATION_LOGIC Cbrt
  #define FUNCTION cbrt
  #include "unaryOverloads.tpp"

  template<typename Real, typename Arg>
  CODI_INLINE PassiveRealType<Real> ceil(ExpressionInterface<Real, Arg> const& arg) {
    return ceil(getPassiveValue(arg.cast()));
  }

  template<typename _Real>
  struct Cos : public UnaryOperation<_Real> {
    public:

      using Real = CODI_DECLARE_DEFAULT(_Real, double);

      /** \copydoc UnaryOperation::primal */
      template<typename Arg>
      static CODI_INLINE Real primal(Arg const& arg) {
        return cos(arg);
      }

      /** \copydoc UnaryOperation::gradient */
      template<typename Arg>
      static CODI_INLINE Real gradient(Arg const& arg, Real const& result) {
        CODI_UNUSED(result);
        return -sin(arg);
      }
  };
  #define OPERATION_LOGIC Cos
  #define FUNCTION cos
  #include "unaryOverloads.tpp"

  template<typename _Real>
  struct Cosh : public UnaryOperation<_Real> {
    public:

      using Real = CODI_DECLARE_DEFAULT(_Real, double);

      /** \copydoc UnaryOperation::primal */
      template<typename Arg>
      static CODI_INLINE Real primal(Arg const& arg) {
        return cosh(arg);
      }

      /** \copydoc UnaryOperation::gradient */
      template<typename Arg>
      static CODI_INLINE Real gradient(Arg const& arg, Real const& result) {
        CODI_UNUSED(result);
        return sinh(arg);
      }
  };
  #define OPERATION_LOGIC Cosh
  #define FUNCTION cosh
  #include "unaryOverloads.tpp"

  template<typename _Real>
  struct Erf : public UnaryOperation<_Real> {
    public:

      using Real = CODI_DECLARE_DEFAULT(_Real, double);

      /** \copydoc UnaryOperation::primal */
      template<typename Arg>
      static CODI_INLINE Real primal(Arg const& arg) {
        return erf(arg);
      }

      /** \copydoc UnaryOperation::gradient */
      template<typename Arg>
      static CODI_INLINE Real gradient(Arg const& arg, Real const& result) {
        CODI_UNUSED(result);
        return 1.128379167095513 * exp( -(arg * arg) ); // erf'(arg) = 2.0 / sqrt(pi) * exp(-arg^2)
      }
  };
  #define OPERATION_LOGIC Erf
  #define FUNCTION erf
  #include "unaryOverloads.tpp"

  template<typename _Real>
  struct Erfc : public UnaryOperation<_Real> {
    public:

      using Real = CODI_DECLARE_DEFAULT(_Real, double);

      /** \copydoc UnaryOperation::primal */
      template<typename Arg>
      static CODI_INLINE Real primal(Arg const& arg) {
        return erfc(arg);
      }

      /** \copydoc UnaryOperation::gradient */
      template<typename Arg>
      static CODI_INLINE Real gradient(Arg const& arg, Real const& result) {
        CODI_UNUSED(result);
        return -1.128379167095513 * exp( -(arg * arg) ); // erfc'(arg) = - 2.0 / sqrt(pi) * exp(-arg^2)
      }
  };
  #define OPERATION_LOGIC Erfc
  #define FUNCTION erfc
  #include "unaryOverloads.tpp"

  template<typename _Real>
  struct Exp : public UnaryOperation<_Real> {
    public:

      using Real = CODI_DECLARE_DEFAULT(_Real, double);

      /** \copydoc UnaryOperation::primal */
      template<typename Arg>
      static CODI_INLINE Real primal(Arg const& arg) {
        return exp(arg);
      }

      /** \copydoc UnaryOperation::gradient */
      template<typename Arg>
      static CODI_INLINE Real gradient(Arg const& arg, Real const& result) {
        CODI_UNUSED(arg);
        return result;
      }
  };
  #define OPERATION_LOGIC Exp
  #define FUNCTION exp
  #include "unaryOverloads.tpp"

  template<typename Real, typename Arg>
  CODI_INLINE PassiveRealType<Real> floor(ExpressionInterface<Real, Arg> const& arg) {
    return floor(getPassiveValue(arg.cast()));
  }

  template<typename Real, typename Arg>
  CODI_INLINE bool isfinite(ExpressionInterface<Real, Arg> const& arg) {
    return isfinite(getPassiveValue(arg.cast()));
  }

  template<typename Real, typename Arg>
  CODI_INLINE bool isinf(ExpressionInterface<Real, Arg> const& arg) {
    return isinf(getPassiveValue(arg.cast()));
  }

  template<typename Real, typename Arg>
  CODI_INLINE bool isnan(ExpressionInterface<Real, Arg> const& arg) {
    return isnan(getPassiveValue(arg.cast()));
  }

  template<typename Real, typename Arg>
  CODI_INLINE bool isnormal(ExpressionInterface<Real, Arg> const& arg) {
    return isnormal(getPassiveValue(arg.cast()));
  }

  template<typename _Real>
  struct Log : public UnaryOperation<_Real> {
    public:

      using Real = CODI_DECLARE_DEFAULT(_Real, double);

      /** \copydoc UnaryOperation::primal */
      template<typename Arg>
      static CODI_INLINE Real primal(Arg const& arg) {
        return log(arg);
      }

      /** \copydoc UnaryOperation::gradient */
      template<typename Arg>
      static CODI_INLINE Real gradient(Arg const& arg, Real const& result) {
        CODI_UNUSED(result);
        if(Config::CheckExpressionArguments) {
          if(0.0 > getPassiveValue(arg)) {
            CODI_EXCEPTION("Logarithm of negative value or zero.(Value: %0.15e)", getPassiveValue(arg));
          }
        }
        return 1.0 / arg;
      }
  };
  #define OPERATION_LOGIC Log
  #define FUNCTION log
  #include "unaryOverloads.tpp"

  template<typename _Real>
  struct Log10 : public UnaryOperation<_Real> {
    public:

      using Real = CODI_DECLARE_DEFAULT(_Real, double);

      /** \copydoc UnaryOperation::primal */
      template<typename Arg>
      static CODI_INLINE Real primal(Arg const& arg) {
        return log10(arg);
      }

      /** \copydoc UnaryOperation::gradient */
      template<typename Arg>
      static CODI_INLINE Real gradient(Arg const& arg, Real const& result) {
        CODI_UNUSED(result);
        if(Config::CheckExpressionArguments) {
          if(0.0 > getPassiveValue(arg)) {
            CODI_EXCEPTION("Logarithm of negative value or zero.(Value: %0.15e)", getPassiveValue(arg));
          }
        }
        return 0.434294481903252 / arg;
      }
  };
  #define OPERATION_LOGIC Log10
  #define FUNCTION log10
  #include "unaryOverloads.tpp"

  template<typename Real, typename Arg>
  CODI_INLINE PassiveRealType<Real> round(ExpressionInterface<Real, Arg> const& arg) {
    using std::round;
    return round(arg.cast().getValue());
  }

  template<typename _Real>
  struct Sin : public UnaryOperation<_Real> {
    public:

      using Real = CODI_DECLARE_DEFAULT(_Real, double);

      /** \copydoc UnaryOperation::primal */
      template<typename Arg>
      static CODI_INLINE Real primal(Arg const& arg) {
        return sin(arg);
      }

      /** \copydoc UnaryOperation::gradient */
      template<typename Arg>
      static CODI_INLINE Real gradient(Arg const& arg, Real const& result) {
        CODI_UNUSED(result);
        return cos(arg);
      }
  };
  #define OPERATION_LOGIC Sin
  #define FUNCTION sin
  #include "unaryOverloads.tpp"

  template<typename _Real>
  struct Sinh : public UnaryOperation<_Real> {
    public:

      using Real = CODI_DECLARE_DEFAULT(_Real, double);

      /** \copydoc UnaryOperation::primal */
      template<typename Arg>
      static CODI_INLINE Real primal(Arg const& arg) {
        return sinh(arg);
      }

      /** \copydoc UnaryOperation::gradient */
      template<typename Arg>
      static CODI_INLINE Real gradient(Arg const& arg, Real const& result) {
        CODI_UNUSED(result);
        return cosh(arg);
      }
  };
  #define OPERATION_LOGIC Sinh
  #define FUNCTION sinh
  #include "unaryOverloads.tpp"

  template<typename _Real>
  struct Sqrt : public UnaryOperation<_Real> {
    public:

      using Real = CODI_DECLARE_DEFAULT(_Real, double);

      /** \copydoc UnaryOperation::primal */
      template<typename Arg>
      static CODI_INLINE Real primal(Arg const& arg) {
        return sqrt(arg);
      }

      /** \copydoc UnaryOperation::gradient */
      template<typename Arg>
      static CODI_INLINE Real gradient(Arg const& arg, Real const& result) {
        if(Config::CheckExpressionArguments) {
          if(0.0 > getPassiveValue(arg)) {
            CODI_EXCEPTION("Sqrt of negative value or zero.(Value: %0.15e)", getPassiveValue(arg));
          }
        }
        if(result != 0.0) {
          return 0.5 / result;
        } else {
          return (Real)0.0;
        }
    }
  };
  #define OPERATION_LOGIC Sqrt
  #define FUNCTION sqrt
  #include "unaryOverloads.tpp"

  template<typename _Real>
  struct Tan : public UnaryOperation<_Real> {
    public:

      using Real = CODI_DECLARE_DEFAULT(_Real, double);

      /** \copydoc UnaryOperation::primal */
      template<typename Arg>
      static CODI_INLINE Real primal(Arg const& arg) {
        return tan(arg);
      }

      /** \copydoc UnaryOperation::gradient */
      template<typename Arg>
      static CODI_INLINE Real gradient(Arg const& arg, Real const& result) {
        CODI_UNUSED(result);
        if(Config::CheckExpressionArguments) {
          if(0.0 == cos(getPassiveValue(arg))) {
            CODI_EXCEPTION("Tan evaluated at (0.5  + i) * PI.(Value: %0.15e)", getPassiveValue(arg));
          }
        }
        Real tmp = 1.0 / cos(arg);
        return tmp * tmp;
      }
  };
  #define OPERATION_LOGIC Tan
  #define FUNCTION tan
  #include "unaryOverloads.tpp"

  template<typename _Real>
  struct Tanh : public UnaryOperation<_Real> {
    public:

      using Real = CODI_DECLARE_DEFAULT(_Real, double);

      /** \copydoc UnaryOperation::primal */
      template<typename Arg>
      static CODI_INLINE Real primal(Arg const& arg) {
        return tanh(arg);
      }

      /** \copydoc UnaryOperation::gradient */
      template<typename Arg>
      static CODI_INLINE Real gradient(Arg const& arg, Real const& result) {
        CODI_UNUSED(arg);
        return 1 - result * result;
      }
  };
  #define OPERATION_LOGIC Tanh
  #define FUNCTION tanh
  #include "unaryOverloads.tpp"

  template<typename _Real>
  struct Tgamma : public UnaryOperation<_Real> {
    public:

      using Real = CODI_DECLARE_DEFAULT(_Real, double);

      /** \copydoc UnaryOperation::primal */
      template<typename Arg>
      static CODI_INLINE Real primal(Arg const& arg) {
        return tgamma(arg);
      }

      /** \copydoc UnaryOperation::gradient */
      template<typename Arg>
      static CODI_INLINE Real gradient(Arg const& arg, Real const& result) {
        if(arg <= 0.0) {
          std::cout << "Derivative for gamma function only for positive arguments at the moment" << std::endl;
          std::exit(1);
        }

        // Implementation of the digamma function is taken from John Burkardt,
        // http://people.sc.fsu.edu/~jburkardt/cpp_src/asa103/asa103.cpp
        //
        // Definition of Gamma(arg): https://en.wikipedia.org/wiki/Gamma_function
        // Definition of DiGamma(arg): https://en.wikipedia.org/wiki/Digamma_function
        // Differentation is Gamma'(arg) = Gamma(arg) * DiGamma(arg)

        Real diGamma = 0.0;
        if(arg <= 0.000001) { // special case for small numbers
          const Real eulerMascheroni = 0.57721566490153286060;
          diGamma = -eulerMascheroni - 1.0/arg + 1.6449340668482264365*arg;
        } else {
          // shift DiGamma(arg) = DiGamma(arg + 1) - 1/arg
          // we require arg large such that the approximation below is more accurate
          Real shiftBound = 8.5;

          Real shiftedValue = arg;
          while( shiftedValue < shiftBound ) {
            diGamma      -= 1.0/shiftedValue;
            shiftedValue += 1.0;
          }

          // Now compute the approximation via an asymptotic series
          Real r = 1.0/shiftedValue;
          diGamma += log(shiftedValue) - 0.5*r;

          Real rSqr = r*r;
          diGamma -= rSqr*(1.0/12.0 - rSqr*(1.0/120.0 - rSqr*(1.0/252.0 - rSqr*(1.0/240.0 - rSqr*(1.0/132.0)))));
        }

        return diGamma*result;
      }
  };
  #define OPERATION_LOGIC Tgamma
  #define FUNCTION tgamma
  #include "unaryOverloads.tpp"
}

namespace std {

  using codi::abs;
  using codi::fabs;
  using codi::acos;
  using codi::asin;
  using codi::atan;
  using codi::atanh;
  using codi::cbrt;
  using codi::ceil;
  using codi::cos;
  using codi::cosh;
  using codi::erf;
  using codi::erfc;
  using codi::exp;
  using codi::floor;
  using codi::isfinite;
  using codi::isinf;
  using codi::isnan;
  using codi::isnormal;
  using codi::log;
  using codi::log10;
  using codi::round;
  using codi::sin;
  using codi::sinh;
  using codi::sqrt;
  using codi::tan;
  using codi::tanh;
  using codi::tgamma;
}
