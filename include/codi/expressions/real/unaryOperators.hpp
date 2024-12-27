/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2025 Chair for Scientific Computing (SciComp), University of Kaiserslautern-Landau
 * Homepage: http://scicomp.rptu.de
 * Contact:  Prof. Nicolas R. Gauger (codi@scicomp.uni-kl.de)
 *
 * Lead developers: Max Sagebaum, Johannes Blühdorn (SciComp, University of Kaiserslautern-Landau)
 *
 * This file is part of CoDiPack (http://scicomp.rptu.de/software/codi).
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

#include <cmath>

#include "../../config.h"
#include "../../misc/exceptions.hpp"
#include "../../misc/macros.hpp"
#include "../../traits/realTraits.hpp"
#include "../unaryExpression.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /*******************************************************************************/
  /// @name Builtin unary operators
  /// @{

  /// UnaryOperation implementation for operator -
  template<typename T_Real>
  struct OperationUnaryMinus : public UnaryOperation<T_Real> {
    public:

      using Real = CODI_DD(T_Real, double);  ///< See BinaryOperation.

      /// \copydoc UnaryOperation::primal
      template<typename Arg>
      static CODI_INLINE Real primal(Arg const& arg) {
        return -arg;
      }

      /// \copydoc UnaryOperation::gradient
      template<typename Arg>
      static CODI_INLINE RealTraits::PassiveReal<Real> gradient(Arg const& arg, Real const& result) {
        CODI_UNUSED(arg);
        CODI_UNUSED(result);
        return -1.0;
      }

      /// \copydoc UnaryOperation::getMathRep()
      static CODI_INLINE std::string getMathRep() {
        return "-";
      }
  };
#define OPERATION_LOGIC OperationUnaryMinus
#define FUNCTION operator-
#include "unaryOverloads.tpp"

  /// Function overload for operator +
  template<typename Real, typename Arg>
  CODI_INLINE ExpressionInterface<Real, Arg> const& operator+(ExpressionInterface<Real, Arg> const& arg) {
    return arg;
  }

  /// @}
  /*******************************************************************************/
  /// @name Standard math library unary operators
  /// @{

  using std::abs;
  using std::acos;
  using std::acosh;
  using std::asin;
  using std::asinh;
  using std::atan;
  using std::atanh;
  using std::cbrt;
  using std::ceil;
  using std::cos;
  using std::cosh;
  using std::erf;
  using std::erfc;
  using std::exp;
  using std::fabs;
  using std::floor;
  using std::isfinite;
  using std::isinf;
  using std::isnan;
  using std::isnormal;
  using std::log;
  using std::log10;
  using std::log1p;
  using std::log2;
  using std::round;
  using std::sin;
  using std::sinh;
  using std::sqrt;
  using std::tan;
  using std::tanh;
  using std::tgamma;

  /// UnaryOperation implementation for abs
  template<typename T_Real>
  struct OperationAbs : public UnaryOperation<T_Real> {
    public:

      using Real = CODI_DD(T_Real, double);  ///< See BinaryOperation.

      /// \copydoc UnaryOperation::primal
      template<typename Arg>
      static CODI_INLINE Real primal(Arg const& arg) {
        return abs(arg);
      }

      /// \copydoc UnaryOperation::gradient
      template<typename Arg>
      static CODI_INLINE RealTraits::PassiveReal<Real> gradient(Arg const& arg, Real const& result) {
        CODI_UNUSED(result);
        if (arg < 0.0) {
          return -1.0;
        } else if (arg > 0.0) {
          return 1.0;
        } else {
          return 0.0;
        }
      }

      /// \copydoc UnaryOperation::getMathRep()
      static CODI_INLINE std::string getMathRep() {
        return "abs";
      }
  };
#define OPERATION_LOGIC OperationAbs
#define FUNCTION abs
#include "unaryOverloads.tpp"

#define OPERATION_LOGIC OperationAbs
#define FUNCTION fabs
#include "unaryOverloads.tpp"

#define OPERATION_LOGIC OperationAbs
#define FUNCTION fabsf
#include "unaryOverloads.tpp"

#define OPERATION_LOGIC OperationAbs
#define FUNCTION fabsl
#include "unaryOverloads.tpp"

  /// UnaryOperation implementation for acos
  template<typename T_Real>
  struct OperationAcos : public UnaryOperation<T_Real> {
    public:

      using Real = CODI_DD(T_Real, double);  ///< See BinaryOperation.

      /// \copydoc UnaryOperation::primal
      template<typename Arg>
      static CODI_INLINE Real primal(Arg const& arg) {
        return acos(arg);
      }

      /// \copydoc UnaryOperation::gradient
      template<typename Arg>
      static CODI_INLINE Real gradient(Arg const& arg, Real const& result) {
        CODI_UNUSED(result);
        if (Config::CheckExpressionArguments) {
          if (RealTraits::getPassiveValue(arg) <= -1.0 || 1.0 <= RealTraits::getPassiveValue(arg)) {
            CODI_EXCEPTION("acos outside of (-1, 1).(Value: %0.15e)", RealTraits::getPassiveValue(arg));
          }
        }
        return -1.0 / sqrt(1.0 - arg * arg);
      }

      /// \copydoc UnaryOperation::getMathRep()
      static CODI_INLINE std::string getMathRep() {
        return "acos";
      }
  };
#define OPERATION_LOGIC OperationAcos
#define FUNCTION acos
#include "unaryOverloads.tpp"

#define OPERATION_LOGIC OperationAcos
#define FUNCTION acosf
#include "unaryOverloads.tpp"

#define OPERATION_LOGIC OperationAcos
#define FUNCTION acosl
#include "unaryOverloads.tpp"

  /// UnaryOperation implementation for acosh
  template<typename T_Real>
  struct OperationAcosh : public UnaryOperation<T_Real> {
    public:

      using Real = CODI_DD(T_Real, double);  ///< See BinaryOperation.

      /// \copydoc UnaryOperation::primal
      template<typename Arg>
      static CODI_INLINE Real primal(Arg const& arg) {
        return acosh(arg);
      }

      /// \copydoc UnaryOperation::gradient
      template<typename Arg>
      static CODI_INLINE Real gradient(Arg const& arg, Real const& result) {
        CODI_UNUSED(result);
        if (Config::CheckExpressionArguments) {
          if (RealTraits::getPassiveValue(arg) <= 1.0) {
            CODI_EXCEPTION("acosh outside of (1, inf).(Value: %0.15e)", RealTraits::getPassiveValue(arg));
          }
        }
        return 1.0 / sqrt(arg * arg - 1.0);
      }

      /// \copydoc UnaryOperation::getMathRep()
      static CODI_INLINE std::string getMathRep() {
        return "acosh";
      }
  };
#define OPERATION_LOGIC OperationAcosh
#define FUNCTION acosh
#include "unaryOverloads.tpp"

#define OPERATION_LOGIC OperationAcosh
#define FUNCTION acoshf
#include "unaryOverloads.tpp"

#define OPERATION_LOGIC OperationAcosh
#define FUNCTION acoshl
#include "unaryOverloads.tpp"

  /// UnaryOperation implementation for asin
  template<typename T_Real>
  struct OperationAsin : public UnaryOperation<T_Real> {
    public:

      using Real = CODI_DD(T_Real, double);  ///< See BinaryOperation.

      /// \copydoc UnaryOperation::primal
      template<typename Arg>
      static CODI_INLINE Real primal(Arg const& arg) {
        return asin(arg);
      }

      /// \copydoc UnaryOperation::gradient
      template<typename Arg>
      static CODI_INLINE Real gradient(Arg const& arg, Real const& result) {
        CODI_UNUSED(result);
        if (Config::CheckExpressionArguments) {
          if (RealTraits::getPassiveValue(arg) <= -1.0 || 1.0 <= RealTraits::getPassiveValue(arg)) {
            CODI_EXCEPTION("asin outside of (-1, 1).(Value: %0.15e)", RealTraits::getPassiveValue(arg));
          }
        }
        return 1.0 / sqrt(1.0 - arg * arg);
      }

      /// \copydoc UnaryOperation::getMathRep()
      static CODI_INLINE std::string getMathRep() {
        return "asin";
      }
  };
#define OPERATION_LOGIC OperationAsin
#define FUNCTION asin
#include "unaryOverloads.tpp"

#define OPERATION_LOGIC OperationAsin
#define FUNCTION asinf
#include "unaryOverloads.tpp"

#define OPERATION_LOGIC OperationAsin
#define FUNCTION asinl
#include "unaryOverloads.tpp"

  /// UnaryOperation implementation for asinh
  template<typename T_Real>
  struct OperationAsinh : public UnaryOperation<T_Real> {
    public:

      using Real = CODI_DD(T_Real, double);  ///< See BinaryOperation.

      /// \copydoc UnaryOperation::primal
      template<typename Arg>
      static CODI_INLINE Real primal(Arg const& arg) {
        return asinh(arg);
      }

      /// \copydoc UnaryOperation::gradient
      template<typename Arg>
      static CODI_INLINE Real gradient(Arg const& arg, Real const& result) {
        CODI_UNUSED(result);
        return 1.0 / sqrt(arg * arg + 1.0);
      }

      /// \copydoc UnaryOperation::getMathRep()
      static CODI_INLINE std::string getMathRep() {
        return "asinh";
      }
  };
#define OPERATION_LOGIC OperationAsinh
#define FUNCTION asinh
#include "unaryOverloads.tpp"

#define OPERATION_LOGIC OperationAsinh
#define FUNCTION asinhf
#include "unaryOverloads.tpp"

#define OPERATION_LOGIC OperationAsinh
#define FUNCTION asinhl
#include "unaryOverloads.tpp"

  /// UnaryOperation implementation for atan
  template<typename T_Real>
  struct OperationAtan : public UnaryOperation<T_Real> {
    public:

      using Real = CODI_DD(T_Real, double);  ///< See BinaryOperation.

      /// \copydoc UnaryOperation::primal
      template<typename Arg>
      static CODI_INLINE Real primal(Arg const& arg) {
        return atan(arg);
      }

      /// \copydoc UnaryOperation::gradient
      template<typename Arg>
      static CODI_INLINE Real gradient(Arg const& arg, Real const& result) {
        CODI_UNUSED(result);
        return 1.0 / (1.0 + arg * arg);
      }

      /// \copydoc UnaryOperation::getMathRep()
      static CODI_INLINE std::string getMathRep() {
        return "atan";
      }
  };
#define OPERATION_LOGIC OperationAtan
#define FUNCTION atan
#include "unaryOverloads.tpp"

#define OPERATION_LOGIC OperationAtan
#define FUNCTION atanf
#include "unaryOverloads.tpp"

#define OPERATION_LOGIC OperationAtan
#define FUNCTION atanl
#include "unaryOverloads.tpp"

  /// UnaryOperation implementation for atanh
  template<typename T_Real>
  struct OperationAtanh : public UnaryOperation<T_Real> {
    public:

      using Real = CODI_DD(T_Real, double);  ///< See BinaryOperation.

      /// \copydoc UnaryOperation::primal
      template<typename Arg>
      static CODI_INLINE Real primal(Arg const& arg) {
        return atanh(arg);
      }

      /// \copydoc UnaryOperation::gradient
      template<typename Arg>
      static CODI_INLINE Real gradient(Arg const& arg, Real const& result) {
        CODI_UNUSED(result);
        if (Config::CheckExpressionArguments) {
          if (RealTraits::getPassiveValue(arg) <= -1.0 || 1.0 <= RealTraits::getPassiveValue(arg)) {
            CODI_EXCEPTION("atanh outside of (-1, 1).(Value: %0.15e)", RealTraits::getPassiveValue(arg));
          }
        }
        return 1.0 / (1.0 - arg * arg);
      }

      /// \copydoc UnaryOperation::getMathRep()
      static CODI_INLINE std::string getMathRep() {
        return "atanh";
      }
  };
#define OPERATION_LOGIC OperationAtanh
#define FUNCTION atanh
#include "unaryOverloads.tpp"

#define OPERATION_LOGIC OperationAtanh
#define FUNCTION atanhf
#include "unaryOverloads.tpp"

#define OPERATION_LOGIC OperationAtanh
#define FUNCTION atanhl
#include "unaryOverloads.tpp"

  /// UnaryOperation implementation for cbrt
  template<typename T_Real>
  struct OperationCbrt : public UnaryOperation<T_Real> {
    public:

      using Real = CODI_DD(T_Real, double);  ///< See BinaryOperation.

      /// \copydoc UnaryOperation::primal
      template<typename Arg>
      static CODI_INLINE Real primal(Arg const& arg) {
        return cbrt(arg);
      }

      /// \copydoc UnaryOperation::gradient
      template<typename Arg>
      static CODI_INLINE Real gradient(Arg const& arg, Real const& result) {
        if (Config::CheckExpressionArguments) {
          if (0.0 == RealTraits::getPassiveValue(arg)) {
            CODI_EXCEPTION("Cbrt of zero value.(Value: %0.15e)", RealTraits::getPassiveValue(arg));
          }
        }
        if (result != 0.0) {
          return 1.0 / (3.0 * result * result);
        } else {
          return (Real)0.0;
        }
      }

      /// \copydoc UnaryOperation::getMathRep()
      static CODI_INLINE std::string getMathRep() {
        return "cbrt";
      }
  };
#define OPERATION_LOGIC OperationCbrt
#define FUNCTION cbrt
#include "unaryOverloads.tpp"

#define OPERATION_LOGIC OperationCbrt
#define FUNCTION cbrtf
#include "unaryOverloads.tpp"

#define OPERATION_LOGIC OperationCbrt
#define FUNCTION cbrtl
#include "unaryOverloads.tpp"

  /// Function overload for ceil
  template<typename Real, typename Arg>
  CODI_INLINE RealTraits::PassiveReal<Real> ceil(ExpressionInterface<Real, Arg> const& arg) {
    return ceil(RealTraits::getPassiveValue(arg.cast()));
  }

  /// Function overload for ceilf
  template<typename Real, typename Arg>
  CODI_INLINE RealTraits::PassiveReal<Real> ceilf(ExpressionInterface<Real, Arg> const& arg) {
    return ceil(RealTraits::getPassiveValue(arg.cast()));
  }

  /// Function overload for ceill
  template<typename Real, typename Arg>
  CODI_INLINE RealTraits::PassiveReal<Real> ceill(ExpressionInterface<Real, Arg> const& arg) {
    return ceil(RealTraits::getPassiveValue(arg.cast()));
  }

  /// UnaryOperation implementation for cos
  template<typename T_Real>
  struct OperationCos : public UnaryOperation<T_Real> {
    public:

      using Real = CODI_DD(T_Real, double);  ///< See BinaryOperation.

      /// \copydoc UnaryOperation::primal
      template<typename Arg>
      static CODI_INLINE Real primal(Arg const& arg) {
        return cos(arg);
      }

      /// \copydoc UnaryOperation::gradient
      template<typename Arg>
      static CODI_INLINE Real gradient(Arg const& arg, Real const& result) {
        CODI_UNUSED(result);
        return -sin(arg);
      }

      /// \copydoc UnaryOperation::getMathRep()
      static CODI_INLINE std::string getMathRep() {
        return "cos";
      }
  };
#define OPERATION_LOGIC OperationCos
#define FUNCTION cos
#include "unaryOverloads.tpp"

#define OPERATION_LOGIC OperationCos
#define FUNCTION cosf
#include "unaryOverloads.tpp"

#define OPERATION_LOGIC OperationCos
#define FUNCTION cosl
#include "unaryOverloads.tpp"

  /// UnaryOperation implementation for cosh
  template<typename T_Real>
  struct OperationCosh : public UnaryOperation<T_Real> {
    public:

      using Real = CODI_DD(T_Real, double);  ///< See BinaryOperation.

      /// \copydoc UnaryOperation::primal
      template<typename Arg>
      static CODI_INLINE Real primal(Arg const& arg) {
        return cosh(arg);
      }

      /// \copydoc UnaryOperation::gradient
      template<typename Arg>
      static CODI_INLINE Real gradient(Arg const& arg, Real const& result) {
        CODI_UNUSED(result);
        return sinh(arg);
      }

      /// \copydoc UnaryOperation::getMathRep()
      static CODI_INLINE std::string getMathRep() {
        return "sinh";
      }
  };
#define OPERATION_LOGIC OperationCosh
#define FUNCTION cosh
#include "unaryOverloads.tpp"

#define OPERATION_LOGIC OperationCosh
#define FUNCTION coshf
#include "unaryOverloads.tpp"

#define OPERATION_LOGIC OperationCosh
#define FUNCTION coshl
#include "unaryOverloads.tpp"

  /// UnaryOperation implementation for erf
  template<typename T_Real>
  struct OperationErf : public UnaryOperation<T_Real> {
    public:

      using Real = CODI_DD(T_Real, double);  ///< See BinaryOperation.

      /// \copydoc UnaryOperation::primal
      template<typename Arg>
      static CODI_INLINE Real primal(Arg const& arg) {
        return erf(arg);
      }

      /// \copydoc UnaryOperation::gradient
      template<typename Arg>
      static CODI_INLINE Real gradient(Arg const& arg, Real const& result) {
        CODI_UNUSED(result);
        return 1.128379167095513 * exp(-(arg * arg));  // erf'(arg) = 2.0 / sqrt(pi) * exp(-arg^2)
      }

      /// \copydoc UnaryOperation::getMathRep()
      static CODI_INLINE std::string getMathRep() {
        return "erf";
      }
  };
#define OPERATION_LOGIC OperationErf
#define FUNCTION erf
#include "unaryOverloads.tpp"

#define OPERATION_LOGIC OperationErf
#define FUNCTION erff
#include "unaryOverloads.tpp"

#define OPERATION_LOGIC OperationErf
#define FUNCTION erfl
#include "unaryOverloads.tpp"

  /// UnaryOperation implementation for erfc
  template<typename T_Real>
  struct OperationErfc : public UnaryOperation<T_Real> {
    public:

      using Real = CODI_DD(T_Real, double);  ///< See BinaryOperation.

      /// \copydoc UnaryOperation::primal
      template<typename Arg>
      static CODI_INLINE Real primal(Arg const& arg) {
        return erfc(arg);
      }

      /// \copydoc UnaryOperation::gradient
      template<typename Arg>
      static CODI_INLINE Real gradient(Arg const& arg, Real const& result) {
        CODI_UNUSED(result);
        return -1.128379167095513 * exp(-(arg * arg));  // erfc'(arg) = - 2.0 / sqrt(pi) * exp(-arg^2)
      }

      /// \copydoc UnaryOperation::getMathRep()
      static CODI_INLINE std::string getMathRep() {
        return "erfc";
      }
  };
#define OPERATION_LOGIC OperationErfc
#define FUNCTION erfc
#include "unaryOverloads.tpp"

#define OPERATION_LOGIC OperationErfc
#define FUNCTION erfcf
#include "unaryOverloads.tpp"

#define OPERATION_LOGIC OperationErfc
#define FUNCTION erfcl
#include "unaryOverloads.tpp"

  /// UnaryOperation implementation for exp
  template<typename T_Real>
  struct OperationExp : public UnaryOperation<T_Real> {
    public:

      using Real = CODI_DD(T_Real, double);  ///< See BinaryOperation.

      /// \copydoc UnaryOperation::primal
      template<typename Arg>
      static CODI_INLINE Real primal(Arg const& arg) {
        return exp(arg);
      }

      /// \copydoc UnaryOperation::gradient
      template<typename Arg>
      static CODI_INLINE Real gradient(Arg const& arg, Real const& result) {
        CODI_UNUSED(arg);
        return result;
      }

      /// \copydoc UnaryOperation::getMathRep()
      static CODI_INLINE std::string getMathRep() {
        return "exp";
      }
  };
#define OPERATION_LOGIC OperationExp
#define FUNCTION exp
#include "unaryOverloads.tpp"

#define OPERATION_LOGIC OperationExp
#define FUNCTION expf
#include "unaryOverloads.tpp"

#define OPERATION_LOGIC OperationExp
#define FUNCTION expl
#include "unaryOverloads.tpp"

  /// Function overload for floor
  template<typename Real, typename Arg>
  CODI_INLINE RealTraits::PassiveReal<Real> floor(ExpressionInterface<Real, Arg> const& arg) {
    return floor(RealTraits::getPassiveValue(arg.cast()));
  }

  /// Function overload for floorf
  template<typename Real, typename Arg>
  CODI_INLINE RealTraits::PassiveReal<Real> floorf(ExpressionInterface<Real, Arg> const& arg) {
    return floor(RealTraits::getPassiveValue(arg.cast()));
  }

  /// Function overload for floorl
  template<typename Real, typename Arg>
  CODI_INLINE RealTraits::PassiveReal<Real> floorl(ExpressionInterface<Real, Arg> const& arg) {
    return floor(RealTraits::getPassiveValue(arg.cast()));
  }

  /// Function overload for isfinite
  template<typename Real, typename Arg>
  CODI_INLINE bool isfinite(ExpressionInterface<Real, Arg> const& arg) {
    return isfinite(RealTraits::getPassiveValue(arg.cast()));
  }

  /// Function overload for isinf
  template<typename Real, typename Arg>
  CODI_INLINE bool isinf(ExpressionInterface<Real, Arg> const& arg) {
    return isinf(RealTraits::getPassiveValue(arg.cast()));
  }

  /// Function overload for isnan
  template<typename Real, typename Arg>
  CODI_INLINE bool isnan(ExpressionInterface<Real, Arg> const& arg) {
    return isnan(RealTraits::getPassiveValue(arg.cast()));
  }

  /// Function overload for isnormal
  template<typename Real, typename Arg>
  CODI_INLINE bool isnormal(ExpressionInterface<Real, Arg> const& arg) {
    return isnormal(RealTraits::getPassiveValue(arg.cast()));
  }

  /// UnaryOperation implementation for log
  template<typename T_Real>
  struct OperationLog : public UnaryOperation<T_Real> {
    public:

      using Real = CODI_DD(T_Real, double);  ///< See BinaryOperation.

      /// \copydoc UnaryOperation::primal
      template<typename Arg>
      static CODI_INLINE Real primal(Arg const& arg) {
        return log(arg);
      }

      /// \copydoc UnaryOperation::gradient
      template<typename Arg>
      static CODI_INLINE Real gradient(Arg const& arg, Real const& result) {
        CODI_UNUSED(result);
        if (Config::CheckExpressionArguments) {
          if (0.0 > RealTraits::getPassiveValue(arg)) {
            CODI_EXCEPTION("Logarithm of negative value or zero.(Value: %0.15e)", RealTraits::getPassiveValue(arg));
          }
        }
        return 1.0 / arg;
      }

      /// \copydoc UnaryOperation::getMathRep()
      static CODI_INLINE std::string getMathRep() {
        return "log";
      }
  };
#define OPERATION_LOGIC OperationLog
#define FUNCTION log
#include "unaryOverloads.tpp"

#define OPERATION_LOGIC OperationLog
#define FUNCTION logf
#include "unaryOverloads.tpp"

#define OPERATION_LOGIC OperationLog
#define FUNCTION logl
#include "unaryOverloads.tpp"

  /// UnaryOperation implementation for log10
  template<typename T_Real>
  struct OperationLog10 : public UnaryOperation<T_Real> {
    public:

      using Real = CODI_DD(T_Real, double);  ///< See BinaryOperation.

      /// \copydoc UnaryOperation::primal
      template<typename Arg>
      static CODI_INLINE Real primal(Arg const& arg) {
        return log10(arg);
      }

      /// \copydoc UnaryOperation::gradient
      template<typename Arg>
      static CODI_INLINE Real gradient(Arg const& arg, Real const& result) {
        CODI_UNUSED(result);
        if (Config::CheckExpressionArguments) {
          if (0.0 > RealTraits::getPassiveValue(arg)) {
            CODI_EXCEPTION("Logarithm of negative value or zero.(Value: %0.15e)", RealTraits::getPassiveValue(arg));
          }
        }
        return 0.434294481903252 / arg;
      }

      /// \copydoc UnaryOperation::getMathRep()
      static CODI_INLINE std::string getMathRep() {
        return "log10";
      }
  };
#define OPERATION_LOGIC OperationLog10
#define FUNCTION log10
#include "unaryOverloads.tpp"

#define OPERATION_LOGIC OperationLog10
#define FUNCTION log10f
#include "unaryOverloads.tpp"

#define OPERATION_LOGIC OperationLog10
#define FUNCTION log10l
#include "unaryOverloads.tpp"

  /// UnaryOperation implementation for log1p
  template<typename T_Real>
  struct OperationLog1p : public UnaryOperation<T_Real> {
    public:

      using Real = CODI_DD(T_Real, double);  ///< See BinaryOperation.

      /// \copydoc UnaryOperation::primal
      template<typename Arg>
      static CODI_INLINE Real primal(Arg const& arg) {
        return log1p(arg);
      }

      /// \copydoc UnaryOperation::gradient
      template<typename Arg>
      static CODI_INLINE Real gradient(Arg const& arg, Real const& result) {
        CODI_UNUSED(result);
        if (Config::CheckExpressionArguments) {
          if (0.0 > RealTraits::getPassiveValue(arg)) {
            CODI_EXCEPTION("Logarithm of negative value or zero.(Value: %0.15e)", RealTraits::getPassiveValue(arg));
          }
        }
        return 1.0 / (arg + 1.0);
      }

      /// \copydoc UnaryOperation::getMathRep()
      static CODI_INLINE std::string getMathRep() {
        return "log1p";
      }
  };
#define OPERATION_LOGIC OperationLog1p
#define FUNCTION log1p
#include "unaryOverloads.tpp"

#define OPERATION_LOGIC OperationLog1p
#define FUNCTION log1pf
#include "unaryOverloads.tpp"

#define OPERATION_LOGIC OperationLog1p
#define FUNCTION log1pl
#include "unaryOverloads.tpp"

  /// UnaryOperation implementation for log2
  template<typename T_Real>
  struct OperationLog2 : public UnaryOperation<T_Real> {
    public:

      using Real = CODI_DD(T_Real, double);  ///< See BinaryOperation.

      /// \copydoc UnaryOperation::primal
      template<typename Arg>
      static CODI_INLINE Real primal(Arg const& arg) {
        return log2(arg);
      }

      /// \copydoc UnaryOperation::gradient
      template<typename Arg>
      static CODI_INLINE Real gradient(Arg const& arg, Real const& result) {
        CODI_UNUSED(result);
        if (Config::CheckExpressionArguments) {
          if (0.0 > RealTraits::getPassiveValue(arg)) {
            CODI_EXCEPTION("Logarithm of negative value or zero.(Value: %0.15e)", RealTraits::getPassiveValue(arg));
          }
        }
        return 1.442695040888963 / arg;
      }

      /// \copydoc UnaryOperation::getMathRep()
      static CODI_INLINE std::string getMathRep() {
        return "log2";
      }
  };
#define OPERATION_LOGIC OperationLog2
#define FUNCTION log2
#include "unaryOverloads.tpp"

#define OPERATION_LOGIC OperationLog2
#define FUNCTION log2f
#include "unaryOverloads.tpp"

#define OPERATION_LOGIC OperationLog2
#define FUNCTION log2l
#include "unaryOverloads.tpp"

  /// Function overload for round
  template<typename Real, typename Arg>
  CODI_INLINE RealTraits::PassiveReal<Real> round(ExpressionInterface<Real, Arg> const& arg) {
    return round(arg.cast().getValue());
  }

  /// Function overload for roundf
  template<typename Real, typename Arg>
  CODI_INLINE RealTraits::PassiveReal<Real> roundf(ExpressionInterface<Real, Arg> const& arg) {
    return round(arg.cast().getValue());
  }

  /// Function overload for roundl
  template<typename Real, typename Arg>
  CODI_INLINE RealTraits::PassiveReal<Real> roundl(ExpressionInterface<Real, Arg> const& arg) {
    return round(arg.cast().getValue());
  }

  /// UnaryOperation implementation for sin
  template<typename T_Real>
  struct OperationSin : public UnaryOperation<T_Real> {
    public:

      using Real = CODI_DD(T_Real, double);  ///< See BinaryOperation.

      /// \copydoc UnaryOperation::primal
      template<typename Arg>
      static CODI_INLINE Real primal(Arg const& arg) {
        return sin(arg);
      }

      /// \copydoc UnaryOperation::gradient
      template<typename Arg>
      static CODI_INLINE Real gradient(Arg const& arg, Real const& result) {
        CODI_UNUSED(result);
        return cos(arg);
      }

      /// \copydoc UnaryOperation::getMathRep()
      static CODI_INLINE std::string getMathRep() {
        return "sin";
      }
  };
#define OPERATION_LOGIC OperationSin
#define FUNCTION sin
#include "unaryOverloads.tpp"

#define OPERATION_LOGIC OperationSin
#define FUNCTION sinf
#include "unaryOverloads.tpp"

#define OPERATION_LOGIC OperationSin
#define FUNCTION sinl
#include "unaryOverloads.tpp"

  /// UnaryOperation implementation for sinh
  template<typename T_Real>
  struct OperationSinh : public UnaryOperation<T_Real> {
    public:

      using Real = CODI_DD(T_Real, double);  ///< See BinaryOperation.

      /// \copydoc UnaryOperation::primal
      template<typename Arg>
      static CODI_INLINE Real primal(Arg const& arg) {
        return sinh(arg);
      }

      /// \copydoc UnaryOperation::gradient
      template<typename Arg>
      static CODI_INLINE Real gradient(Arg const& arg, Real const& result) {
        CODI_UNUSED(result);
        return cosh(arg);
      }

      /// \copydoc UnaryOperation::getMathRep()
      static CODI_INLINE std::string getMathRep() {
        return "sinh";
      }
  };
#define OPERATION_LOGIC OperationSinh
#define FUNCTION sinh
#include "unaryOverloads.tpp"

#define OPERATION_LOGIC OperationSinh
#define FUNCTION sinhf
#include "unaryOverloads.tpp"

#define OPERATION_LOGIC OperationSinh
#define FUNCTION sinhl
#include "unaryOverloads.tpp"

  /// UnaryOperation implementation for sqrt
  template<typename T_Real>
  struct OperationSqrt : public UnaryOperation<T_Real> {
    public:

      using Real = CODI_DD(T_Real, double);  ///< See BinaryOperation.

      /// \copydoc UnaryOperation::primal
      template<typename Arg>
      static CODI_INLINE Real primal(Arg const& arg) {
        return sqrt(arg);
      }

      /// \copydoc UnaryOperation::gradient
      template<typename Arg>
      static CODI_INLINE Real gradient(Arg const& arg, Real const& result) {
        if (Config::CheckExpressionArguments) {
          if (0.0 > RealTraits::getPassiveValue(arg)) {
            CODI_EXCEPTION("Sqrt of negative value or zero.(Value: %0.15e)", RealTraits::getPassiveValue(arg));
          }
        }
        if (result != 0.0) {
          return 0.5 / result;
        } else {
          return (Real)0.0;
        }
      }

      /// \copydoc UnaryOperation::getMathRep()
      static CODI_INLINE std::string getMathRep() {
        return "sqrt";
      }
  };
#define OPERATION_LOGIC OperationSqrt
#define FUNCTION sqrt
#include "unaryOverloads.tpp"

#define OPERATION_LOGIC OperationSqrt
#define FUNCTION sqrtf
#include "unaryOverloads.tpp"

#define OPERATION_LOGIC OperationSqrt
#define FUNCTION sqrtl
#include "unaryOverloads.tpp"

  /// UnaryOperation implementation for tan
  template<typename T_Real>
  struct OperationTan : public UnaryOperation<T_Real> {
    public:

      using Real = CODI_DD(T_Real, double);  ///< See BinaryOperation.

      /// \copydoc UnaryOperation::primal
      template<typename Arg>
      static CODI_INLINE Real primal(Arg const& arg) {
        return tan(arg);
      }

      /// \copydoc UnaryOperation::gradient
      template<typename Arg>
      static CODI_INLINE Real gradient(Arg const& arg, Real const& result) {
        CODI_UNUSED(result);
        if (Config::CheckExpressionArguments) {
          if (0.0 == cos(RealTraits::getPassiveValue(arg))) {
            CODI_EXCEPTION("Tan evaluated at (0.5  + i) * PI.(Value: %0.15e)", RealTraits::getPassiveValue(arg));
          }
        }
        Real tmp = 1.0 / cos(arg);
        return tmp * tmp;
      }

      /// \copydoc UnaryOperation::getMathRep()
      static CODI_INLINE std::string getMathRep() {
        return "tan";
      }
  };
#define OPERATION_LOGIC OperationTan
#define FUNCTION tan
#include "unaryOverloads.tpp"

#define OPERATION_LOGIC OperationTan
#define FUNCTION tanf
#include "unaryOverloads.tpp"

#define OPERATION_LOGIC OperationTan
#define FUNCTION tanl
#include "unaryOverloads.tpp"

  /// UnaryOperation implementation for tanh
  template<typename T_Real>
  struct OperationTanh : public UnaryOperation<T_Real> {
    public:

      using Real = CODI_DD(T_Real, double);  ///< See BinaryOperation.

      /// \copydoc UnaryOperation::primal
      template<typename Arg>
      static CODI_INLINE Real primal(Arg const& arg) {
        return tanh(arg);
      }

      /// \copydoc UnaryOperation::gradient
      template<typename Arg>
      static CODI_INLINE Real gradient(Arg const& arg, Real const& result) {
        CODI_UNUSED(arg);
        return 1.0 - result * result;
      }

      /// \copydoc UnaryOperation::getMathRep()
      static CODI_INLINE std::string getMathRep() {
        return "tanh";
      }
  };
#define OPERATION_LOGIC OperationTanh
#define FUNCTION tanh
#include "unaryOverloads.tpp"

#define OPERATION_LOGIC OperationTanh
#define FUNCTION tanhf
#include "unaryOverloads.tpp"

#define OPERATION_LOGIC OperationTanh
#define FUNCTION tanhl
#include "unaryOverloads.tpp"

  /// UnaryOperation implementation for tgamma
  template<typename T_Real>
  struct OperationTgamma : public UnaryOperation<T_Real> {
    public:

      using Real = CODI_DD(T_Real, double);  ///< See BinaryOperation.

      /// \copydoc UnaryOperation::primal
      template<typename Arg>
      static CODI_INLINE Real primal(Arg const& arg) {
        return tgamma(arg);
      }

      /// \copydoc UnaryOperation::gradient
      template<typename Arg>
      static CODI_INLINE Real gradient(Arg const& arg, Real const& result) {
        if (arg <= 0.0) {
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
        if (arg <= 0.000001) {  // Special case for small numbers.
          const Real eulerMascheroni = 0.57721566490153286060;
          diGamma = -eulerMascheroni - 1.0 / arg + 1.6449340668482264365 * arg;
        } else {
          // Shift DiGamma(arg) = DiGamma(arg + 1) - 1/arg
          // We require arg large such that the approximation below is more accurate.
          Real shiftBound = 8.5;

          Real shiftedValue = arg;
          while (shiftedValue < shiftBound) {
            diGamma -= 1.0 / shiftedValue;
            shiftedValue += 1.0;
          }

          // Now compute the approximation via an asymptotic series.
          Real r = 1.0 / shiftedValue;
          diGamma += log(shiftedValue) - 0.5 * r;

          Real rSqr = r * r;
          diGamma -= rSqr * (1.0 / 12.0 -
                             rSqr * (1.0 / 120.0 - rSqr * (1.0 / 252.0 - rSqr * (1.0 / 240.0 - rSqr * (1.0 / 132.0)))));
        }

        return diGamma * result;
      }

      /// \copydoc UnaryOperation::getMathRep()
      static CODI_INLINE std::string getMathRep() {
        return "tgamma";
      }
  };
#define OPERATION_LOGIC OperationTgamma
#define FUNCTION tgamma
#include "unaryOverloads.tpp"

#define OPERATION_LOGIC OperationTgamma
#define FUNCTION tgammaf
#include "unaryOverloads.tpp"

#define OPERATION_LOGIC OperationTgamma
#define FUNCTION tgammal
#include "unaryOverloads.tpp"

  /// @}
  /*******************************************************************************/
  /// @name Additional standard library unary operators
  /// @{

  /// Function overload for to_string
  template<typename Real, typename Arg>
  CODI_INLINE std::string to_string(ExpressionInterface<Real, Arg> const& arg) {
    using std::to_string;

    return to_string(RealTraits::getPassiveValue(arg.cast()));
  }

  /// @}
}

namespace std {

  using codi::abs;
  using codi::acos;
  using codi::acosf;
  using codi::acosh;
  using codi::acoshf;
  using codi::acoshl;
  using codi::acosl;
  using codi::asin;
  using codi::asinf;
  using codi::asinh;
  using codi::asinhf;
  using codi::asinhl;
  using codi::asinl;
  using codi::atan;
  using codi::atanf;
  using codi::atanh;
  using codi::atanhf;
  using codi::atanhl;
  using codi::atanl;
  using codi::cbrt;
  using codi::cbrtf;
  using codi::cbrtl;
  using codi::ceil;
  using codi::ceilf;
  using codi::ceill;
  using codi::cos;
  using codi::cosf;
  using codi::cosh;
  using codi::coshf;
  using codi::coshl;
  using codi::cosl;
  using codi::erf;
  using codi::erfc;
  using codi::erfcf;
  using codi::erfcl;
  using codi::erff;
  using codi::erfl;
  using codi::exp;
  using codi::expf;
  using codi::expl;
  using codi::fabs;
  using codi::fabsf;
  using codi::fabsl;
  using codi::floor;
  using codi::floorf;
  using codi::floorl;
  using codi::isfinite;
  using codi::isinf;
  using codi::isnan;
  using codi::isnormal;
  using codi::log;
  using codi::log10;
  using codi::log10f;
  using codi::log10l;
  using codi::log1p;
  using codi::log1pf;
  using codi::log1pl;
  using codi::log2;
  using codi::log2f;
  using codi::log2l;
  using codi::logf;
  using codi::logl;
  using codi::round;
  using codi::roundf;
  using codi::roundl;
  using codi::sin;
  using codi::sinf;
  using codi::sinh;
  using codi::sinhf;
  using codi::sinhl;
  using codi::sinl;
  using codi::sqrt;
  using codi::sqrtf;
  using codi::sqrtl;
  using codi::tan;
  using codi::tanf;
  using codi::tanh;
  using codi::tanhf;
  using codi::tanhl;
  using codi::tanl;
  using codi::tgamma;
  using codi::tgammaf;
  using codi::tgammal;

  using codi::to_string;
}
