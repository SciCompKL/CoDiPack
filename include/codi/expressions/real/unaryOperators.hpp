/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2026 Chair for Scientific Computing (SciComp), RPTU University Kaiserslautern-Landau
 * Homepage: http://scicomp.rptu.de
 * Contact:  Prof. Nicolas R. Gauger (codi@scicomp.uni-kl.de)
 *
 * Lead developers: Max Sagebaum, Johannes Blühdorn (SciComp, RPTU University Kaiserslautern-Landau)
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
 *  - SciComp, RPTU University Kaiserslautern-Landau:
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
#include "../computeExpression.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /*******************************************************************************/
  /// @name Builtin unary operators
  /// @{

  /// UnaryJacobianOperation implementation for operator -
  template<typename T_Real>
  struct OperationUnaryMinus : public UnaryJacobianOperation<T_Real, OperationUnaryMinus<T_Real>> {
    public:

      using Real = CODI_DD(T_Real, double);  ///< See UnaryJacobianOperation.

      /// \copydoc UnaryJacobianOperation::primal
      template<typename Arg>
      static CODI_INLINE Real primal(Arg const& arg) {
        return -arg;
      }

      /// \copydoc UnaryJacobianOperation::gradient
      template<typename Arg>
      static CODI_INLINE RealTraits::PassiveReal<Real> gradient(Arg const& arg, Real const& result) {
        CODI_UNUSED(arg);
        CODI_UNUSED(result);
        return -1.0;
      }

      /// \copydoc UnaryJacobianOperation::getMathRep()
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

  /// UnaryJacobianOperation implementation for abs
  template<typename T_Real>
  struct OperationAbs : public UnaryJacobianOperation<T_Real, OperationAbs<T_Real>> {
    public:

      using Real = CODI_DD(T_Real, double);  ///< See UnaryJacobianOperation.

      /// \copydoc UnaryJacobianOperation::primal
      template<typename Arg>
      static CODI_INLINE Real primal(Arg const& arg) {
        return abs(arg);
      }

      /// \copydoc UnaryJacobianOperation::gradient
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

      /// \copydoc UnaryJacobianOperation::getMathRep()
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

  /// UnaryJacobianOperation implementation for acos
  template<typename T_Real>
  struct OperationAcos : public UnaryJacobianOperation<T_Real, OperationAcos<T_Real>> {
    public:

      using Real = CODI_DD(T_Real, double);  ///< See UnaryJacobianOperation.

      /// \copydoc UnaryJacobianOperation::primal
      template<typename Arg>
      static CODI_INLINE Real primal(Arg const& arg) {
        return acos(arg);
      }

      /// \copydoc UnaryJacobianOperation::gradient
      template<typename Arg>
      static CODI_INLINE Real gradient(Arg const& arg, Real const& result) {
        CODI_UNUSED(result);
        if (Config::CheckExpressionArguments) {
          checkArgument(arg);
        }
        return -1.0 / sqrt(1.0 - arg * arg);
      }

      /// \copydoc UnaryJacobianOperation::getMathRep()
      static CODI_INLINE std::string getMathRep() {
        return "acos";
      }

    private:

      template<typename Arg>
      CODI_INLINE static void checkArgument(Arg const& arg) {
        if (RealTraits::getPassiveValue(arg) <= -1.0 || 1.0 <= RealTraits::getPassiveValue(arg)) {
          CODI_EXCEPTION("Argument of acos outside of (-1, 1). (Value: %0.15e)", RealTraits::getPassiveValue(arg));
        }
      }

      template<typename Arg>
      CODI_INLINE static void checkArgument(std::complex<Arg> const& arg) {
        if (0.0 == RealTraits::getPassiveValue(imag(arg)) &&
            (RealTraits::getPassiveValue(real(arg)) == -1.0 || 1.0 == RealTraits::getPassiveValue(real(arg)))) {
          CODI_EXCEPTION("Argument of acos outside of C \\ {-1, 1}. (Value: %0.15e + %0.15e i)",
                         RealTraits::getPassiveValue(real(arg)), RealTraits::getPassiveValue(imag(arg)));
        }
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

  /// UnaryJacobianOperation implementation for acosh
  template<typename T_Real>
  struct OperationAcosh : public UnaryJacobianOperation<T_Real, OperationAcosh<T_Real>> {
    public:

      using Real = CODI_DD(T_Real, double);  ///< See UnaryJacobianOperation.

      /// \copydoc UnaryJacobianOperation::primal
      template<typename Arg>
      static CODI_INLINE Real primal(Arg const& arg) {
        return acosh(arg);
      }

      /// \copydoc UnaryJacobianOperation::gradient
      template<typename Arg>
      static CODI_INLINE Real gradient(Arg const& arg, Real const& result) {
        CODI_UNUSED(result);
        if (Config::CheckExpressionArguments) {
          checkArgument(arg);
        }
        return 1.0 / (sqrt(arg + 1.0) * sqrt(arg - 1.0));
      }

      /// \copydoc UnaryJacobianOperation::getMathRep()
      static CODI_INLINE std::string getMathRep() {
        return "acosh";
      }

    private:

      template<typename Arg>
      CODI_INLINE static void checkArgument(Arg const& arg) {
        if (RealTraits::getPassiveValue(arg) > 1.0) {
          CODI_EXCEPTION("Argument of acosh outside of (1, inf). (Value: %0.15e)", RealTraits::getPassiveValue(arg));
        }
      }

      template<typename Arg>
      CODI_INLINE static void checkArgument(std::complex<Arg> const& arg) {
        if (0.0 == RealTraits::getPassiveValue(imag(arg)) && RealTraits::getPassiveValue(real(arg)) >= 1.0) {
          CODI_EXCEPTION("Argument of acosh outside of C \\ {1, -1}. (Value: %0.15e + %0.15e i)",
                         RealTraits::getPassiveValue(real(arg)), RealTraits::getPassiveValue(imag(arg)));
        }
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

  /// UnaryJacobianOperation implementation for asin
  template<typename T_Real>
  struct OperationAsin : public UnaryJacobianOperation<T_Real, OperationAsin<T_Real>> {
    public:

      using Real = CODI_DD(T_Real, double);  ///< See UnaryJacobianOperation.

      /// \copydoc UnaryJacobianOperation::primal
      template<typename Arg>
      static CODI_INLINE Real primal(Arg const& arg) {
        return asin(arg);
      }

      /// \copydoc UnaryJacobianOperation::gradient
      template<typename Arg>
      static CODI_INLINE Real gradient(Arg const& arg, Real const& result) {
        CODI_UNUSED(result);
        if (Config::CheckExpressionArguments) {
          checkArgument(arg);
        }
        return 1.0 / sqrt(1.0 - arg * arg);
      }

      /// \copydoc UnaryJacobianOperation::getMathRep()
      static CODI_INLINE std::string getMathRep() {
        return "asin";
      }

    private:

      template<typename Arg>
      CODI_INLINE static void checkArgument(Arg const& arg) {
        if (RealTraits::getPassiveValue(arg) <= -1.0 || 1.0 <= RealTraits::getPassiveValue(arg)) {
          CODI_EXCEPTION("Argument of asin outside of (-1, 1). (Value: %0.15e)", RealTraits::getPassiveValue(arg));
        }
      }

      template<typename Arg>
      CODI_INLINE static void checkArgument(std::complex<Arg> const& arg) {
        if (0.0 == RealTraits::getPassiveValue(imag(arg)) &&
            (RealTraits::getPassiveValue(real(arg)) == -1.0 || 1.0 == RealTraits::getPassiveValue(real(arg)))) {
          CODI_EXCEPTION("Argument of asin outside of C \\ {1, -1}. (Value: %0.15e + %0.15e i)",
                         RealTraits::getPassiveValue(real(arg)), RealTraits::getPassiveValue(imag(arg)));
        }
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

  /// UnaryJacobianOperation implementation for asinh
  template<typename T_Real>
  struct OperationAsinh : public UnaryJacobianOperation<T_Real, OperationAsinh<T_Real>> {
    public:

      using Real = CODI_DD(T_Real, double);  ///< See UnaryJacobianOperation.

      /// \copydoc UnaryJacobianOperation::primal
      template<typename Arg>
      static CODI_INLINE Real primal(Arg const& arg) {
        return asinh(arg);
      }

      /// \copydoc UnaryJacobianOperation::gradient
      template<typename Arg>
      static CODI_INLINE Real gradient(Arg const& arg, Real const& result) {
        CODI_UNUSED(result);

        if (Config::CheckExpressionArguments) {
          checkArgument(arg);
        }

        return 1.0 / sqrt(arg * arg + 1.0);
      }

      /// \copydoc UnaryJacobianOperation::getMathRep()
      static CODI_INLINE std::string getMathRep() {
        return "asinh";
      }

    private:

      template<typename Arg>
      CODI_INLINE static void checkArgument(Arg const& arg) {
        CODI_UNUSED(arg);
        // Nothing to check.
      }

      template<typename Arg>
      CODI_INLINE static void checkArgument(std::complex<Arg> const& arg) {
        if (1.0 == abs(RealTraits::getPassiveValue(imag(arg))) && (0 == RealTraits::getPassiveValue(real(arg)))) {
          CODI_EXCEPTION("Argument of asinh outside of C \\ {i, -i}. (Value: %0.15e + %0.15e i)",
                         RealTraits::getPassiveValue(real(arg)), RealTraits::getPassiveValue(imag(arg)));
        }
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

  /// UnaryJacobianOperation implementation for atan
  template<typename T_Real>
  struct OperationAtan : public UnaryJacobianOperation<T_Real, OperationAtan<T_Real>> {
    public:

      using Real = CODI_DD(T_Real, double);  ///< See UnaryJacobianOperation.

      /// \copydoc UnaryJacobianOperation::primal
      template<typename Arg>
      static CODI_INLINE Real primal(Arg const& arg) {
        return atan(arg);
      }

      /// \copydoc UnaryJacobianOperation::gradient
      template<typename Arg>
      static CODI_INLINE Real gradient(Arg const& arg, Real const& result) {
        CODI_UNUSED(result);

        if (Config::CheckExpressionArguments) {
          checkArgument(arg);
        }

        return 1.0 / (1.0 + arg * arg);
      }

      /// \copydoc UnaryJacobianOperation::getMathRep()
      static CODI_INLINE std::string getMathRep() {
        return "atan";
      }

    private:

      template<typename Arg>
      CODI_INLINE static void checkArgument(Arg const& arg) {
        CODI_UNUSED(arg);
        // Nothing to check.
      }

      template<typename Arg>
      CODI_INLINE static void checkArgument(std::complex<Arg> const& arg) {
        if (1.0 == abs(RealTraits::getPassiveValue(imag(arg))) && (0 == RealTraits::getPassiveValue(real(arg)))) {
          CODI_EXCEPTION("Argument of atan outside of C \\ {i, -i}. (Value: %0.15e + %0.15e i)",
                         RealTraits::getPassiveValue(real(arg)), RealTraits::getPassiveValue(imag(arg)));
        }
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

  /// UnaryJacobianOperation implementation for atanh
  template<typename T_Real>
  struct OperationAtanh : public UnaryJacobianOperation<T_Real, OperationAtanh<T_Real>> {
    public:

      using Real = CODI_DD(T_Real, double);  ///< See UnaryJacobianOperation.

      /// \copydoc UnaryJacobianOperation::primal
      template<typename Arg>
      static CODI_INLINE Real primal(Arg const& arg) {
        return atanh(arg);
      }

      /// \copydoc UnaryJacobianOperation::gradient
      template<typename Arg>
      static CODI_INLINE Real gradient(Arg const& arg, Real const& result) {
        CODI_UNUSED(result);
        if (Config::CheckExpressionArguments) {
          checkArgument(arg);
        }
        return 1.0 / (1.0 - arg * arg);
      }

      /// \copydoc UnaryJacobianOperation::getMathRep()
      static CODI_INLINE std::string getMathRep() {
        return "atanh";
      }

    private:

      template<typename Arg>
      CODI_INLINE static void checkArgument(Arg const& arg) {
        if (RealTraits::getPassiveValue(arg) <= -1.0 || 1.0 <= RealTraits::getPassiveValue(arg)) {
          CODI_EXCEPTION("Argument of atanh outside of (-1, 1). (Value: %0.15e)", RealTraits::getPassiveValue(arg));
        }
      }

      template<typename Arg>
      CODI_INLINE static void checkArgument(std::complex<Arg> const& arg) {
        if (0.0 == RealTraits::getPassiveValue(imag(arg)) &&
            (RealTraits::getPassiveValue(real(arg)) <= -1.0 || 1.0 <= RealTraits::getPassiveValue(real(arg)))) {
          CODI_EXCEPTION("Argument of atanh outside of C \\ (-1, 1). (Value: %0.15e + %0.15e i)",
                         RealTraits::getPassiveValue(real(arg)), RealTraits::getPassiveValue(imag(arg)));
        }
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

  /// UnaryJacobianOperation implementation for cbrt
  template<typename T_Real>
  struct OperationCbrt : public UnaryJacobianOperation<T_Real, OperationCbrt<T_Real>> {
    public:

      using Real = CODI_DD(T_Real, double);  ///< See UnaryJacobianOperation.

      /// \copydoc UnaryJacobianOperation::primal
      template<typename Arg>
      static CODI_INLINE Real primal(Arg const& arg) {
        return cbrt(arg);
      }

      /// \copydoc UnaryJacobianOperation::gradient
      template<typename Arg>
      static CODI_INLINE Real gradient(Arg const& arg, Real const& result) {
        if (Config::CheckExpressionArguments) {
          if (0.0 == RealTraits::getPassiveValue(arg)) {
            CODI_EXCEPTION("Cbrt of zero value. (Value: %0.15e)", RealTraits::getPassiveValue(arg));
          }
        }
        if (result != 0.0) {
          return 1.0 / (3.0 * result * result);
        } else {
          return (Real)0.0;
        }
      }

      /// \copydoc UnaryJacobianOperation::getMathRep()
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

  /// UnaryJacobianOperation implementation for cos
  template<typename T_Real>
  struct OperationCos : public UnaryJacobianOperation<T_Real, OperationCos<T_Real>> {
    public:

      using Real = CODI_DD(T_Real, double);  ///< See UnaryJacobianOperation.

      /// \copydoc UnaryJacobianOperation::primal
      template<typename Arg>
      static CODI_INLINE Real primal(Arg const& arg) {
        return cos(arg);
      }

      /// \copydoc UnaryJacobianOperation::gradient
      template<typename Arg>
      static CODI_INLINE Real gradient(Arg const& arg, Real const& result) {
        CODI_UNUSED(result);
        return -sin(arg);
      }

      /// \copydoc UnaryJacobianOperation::getMathRep()
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

  /// UnaryJacobianOperation implementation for cosh
  template<typename T_Real>
  struct OperationCosh : public UnaryJacobianOperation<T_Real, OperationCosh<T_Real>> {
    public:

      using Real = CODI_DD(T_Real, double);  ///< See UnaryJacobianOperation.

      /// \copydoc UnaryJacobianOperation::primal
      template<typename Arg>
      static CODI_INLINE Real primal(Arg const& arg) {
        return cosh(arg);
      }

      /// \copydoc UnaryJacobianOperation::gradient
      template<typename Arg>
      static CODI_INLINE Real gradient(Arg const& arg, Real const& result) {
        CODI_UNUSED(result);
        return sinh(arg);
      }

      /// \copydoc UnaryJacobianOperation::getMathRep()
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

  /// UnaryJacobianOperation implementation for erf
  template<typename T_Real>
  struct OperationErf : public UnaryJacobianOperation<T_Real, OperationErf<T_Real>> {
    public:

      using Real = CODI_DD(T_Real, double);  ///< See UnaryJacobianOperation.

      /// \copydoc UnaryJacobianOperation::primal
      template<typename Arg>
      static CODI_INLINE Real primal(Arg const& arg) {
        return erf(arg);
      }

      /// \copydoc UnaryJacobianOperation::gradient
      template<typename Arg>
      static CODI_INLINE Real gradient(Arg const& arg, Real const& result) {
        CODI_UNUSED(result);
        return 1.128379167095513 * exp(-(arg * arg));  // erf'(arg) = 2.0 / sqrt(pi) * exp(-arg^2)
      }

      /// \copydoc UnaryJacobianOperation::getMathRep()
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

  /// UnaryJacobianOperation implementation for erfc
  template<typename T_Real>
  struct OperationErfc : public UnaryJacobianOperation<T_Real, OperationErfc<T_Real>> {
    public:

      using Real = CODI_DD(T_Real, double);  ///< See UnaryJacobianOperation.

      /// \copydoc UnaryJacobianOperation::primal
      template<typename Arg>
      static CODI_INLINE Real primal(Arg const& arg) {
        return erfc(arg);
      }

      /// \copydoc UnaryJacobianOperation::gradient
      template<typename Arg>
      static CODI_INLINE Real gradient(Arg const& arg, Real const& result) {
        CODI_UNUSED(result);
        return -1.128379167095513 * exp(-(arg * arg));  // erfc'(arg) = - 2.0 / sqrt(pi) * exp(-arg^2)
      }

      /// \copydoc UnaryJacobianOperation::getMathRep()
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

  /// UnaryJacobianOperation implementation for exp
  template<typename T_Real>
  struct OperationExp : public UnaryJacobianOperation<T_Real, OperationExp<T_Real>> {
    public:

      using Real = CODI_DD(T_Real, double);  ///< See UnaryJacobianOperation.

      /// \copydoc UnaryJacobianOperation::primal
      template<typename Arg>
      static CODI_INLINE Real primal(Arg const& arg) {
        return exp(arg);
      }

      /// \copydoc UnaryJacobianOperation::gradient
      template<typename Arg>
      static CODI_INLINE Real gradient(Arg const& arg, Real const& result) {
        CODI_UNUSED(arg);
        return result;
      }

      /// \copydoc UnaryJacobianOperation::getMathRep()
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

  /// UnaryJacobianOperation implementation for log
  template<typename T_Real>
  struct OperationLog : public UnaryJacobianOperation<T_Real, OperationLog<T_Real>> {
    public:

      using Real = CODI_DD(T_Real, double);  ///< See UnaryJacobianOperation.

      /// \copydoc UnaryJacobianOperation::primal
      template<typename Arg>
      static CODI_INLINE Real primal(Arg const& arg) {
        return log(arg);
      }

      /// \copydoc UnaryJacobianOperation::gradient
      template<typename Arg>
      static CODI_INLINE Real gradient(Arg const& arg, Real const& result) {
        CODI_UNUSED(result);
        if (Config::CheckExpressionArguments) {
          checkArgument(arg);
        }
        return 1.0 / arg;
      }

      /// \copydoc UnaryJacobianOperation::getMathRep()
      static CODI_INLINE std::string getMathRep() {
        return "log";
      }

    private:

      template<typename Arg>
      CODI_INLINE static void checkArgument(Arg const& arg) {
        if (0.0 > RealTraits::getPassiveValue(arg)) {
          CODI_EXCEPTION("Logarithm of negative value or zero. (Value: %0.15e)", RealTraits::getPassiveValue(arg));
        }
      }

      template<typename Arg>
      CODI_INLINE static void checkArgument(std::complex<Arg> const& arg) {
        if (0.0 == RealTraits::getPassiveValue(abs(arg))) {
          CODI_EXCEPTION("Logarithm of zero. (Value: %0.15e)", RealTraits::getPassiveValue(abs(arg)));
        }
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

  /// UnaryJacobianOperation implementation for log10
  template<typename T_Real>
  struct OperationLog10 : public UnaryJacobianOperation<T_Real, OperationLog10<T_Real>> {
    public:

      using Real = CODI_DD(T_Real, double);  ///< See UnaryJacobianOperation.

      /// \copydoc UnaryJacobianOperation::primal
      template<typename Arg>
      static CODI_INLINE Real primal(Arg const& arg) {
        return log10(arg);
      }

      /// \copydoc UnaryJacobianOperation::gradient
      template<typename Arg>
      static CODI_INLINE Real gradient(Arg const& arg, Real const& result) {
        CODI_UNUSED(result);
        if (Config::CheckExpressionArguments) {
          checkArgument(arg);
        }
        return 0.434294481903252 / arg;
      }

      /// \copydoc UnaryJacobianOperation::getMathRep()
      static CODI_INLINE std::string getMathRep() {
        return "log10";
      }

    private:
      template<typename Arg>
      CODI_INLINE static void checkArgument(Arg const& arg) {
        if (0.0 > RealTraits::getPassiveValue(arg)) {
          CODI_EXCEPTION("Logarithm of negative value or zero. (Value: %0.15e)", RealTraits::getPassiveValue(arg));
        }
      }

      template<typename Arg>
      CODI_INLINE static void checkArgument(std::complex<Arg> const& arg) {
        if (0.0 == RealTraits::getPassiveValue(abs(arg))) {
          CODI_EXCEPTION("Logarithm of zero. (Value: %0.15e)", RealTraits::getPassiveValue(abs(arg)));
        }
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

  /// UnaryJacobianOperation implementation for log1p
  template<typename T_Real>
  struct OperationLog1p : public UnaryJacobianOperation<T_Real, OperationLog1p<T_Real>> {
    public:

      using Real = CODI_DD(T_Real, double);  ///< See UnaryJacobianOperation.

      /// \copydoc UnaryJacobianOperation::primal
      template<typename Arg>
      static CODI_INLINE Real primal(Arg const& arg) {
        return log1p(arg);
      }

      /// \copydoc UnaryJacobianOperation::gradient
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

      /// \copydoc UnaryJacobianOperation::getMathRep()
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

  /// UnaryJacobianOperation implementation for log2
  template<typename T_Real>
  struct OperationLog2 : public UnaryJacobianOperation<T_Real, OperationLog2<T_Real>> {
    public:

      using Real = CODI_DD(T_Real, double);  ///< See UnaryJacobianOperation.

      /// \copydoc UnaryJacobianOperation::primal
      template<typename Arg>
      static CODI_INLINE Real primal(Arg const& arg) {
        return log2(arg);
      }

      /// \copydoc UnaryJacobianOperation::gradient
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

      /// \copydoc UnaryJacobianOperation::getMathRep()
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

  /// UnaryJacobianOperation implementation for sin
  template<typename T_Real>
  struct OperationSin : public UnaryJacobianOperation<T_Real, OperationSin<T_Real>> {
    public:

      using Real = CODI_DD(T_Real, double);  ///< See UnaryJacobianOperation.

      /// \copydoc UnaryJacobianOperation::primal
      template<typename Arg>
      static CODI_INLINE Real primal(Arg const& arg) {
        return sin(arg);
      }

      /// \copydoc UnaryJacobianOperation::gradient
      template<typename Arg>
      static CODI_INLINE Real gradient(Arg const& arg, Real const& result) {
        CODI_UNUSED(result);
        return cos(arg);
      }

      /// \copydoc UnaryJacobianOperation::getMathRep()
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

  /// UnaryJacobianOperation implementation for sinh
  template<typename T_Real>
  struct OperationSinh : public UnaryJacobianOperation<T_Real, OperationSinh<T_Real>> {
    public:

      using Real = CODI_DD(T_Real, double);  ///< See UnaryJacobianOperation.

      /// \copydoc UnaryJacobianOperation::primal
      template<typename Arg>
      static CODI_INLINE Real primal(Arg const& arg) {
        return sinh(arg);
      }

      /// \copydoc UnaryJacobianOperation::gradient
      template<typename Arg>
      static CODI_INLINE Real gradient(Arg const& arg, Real const& result) {
        CODI_UNUSED(result);
        return cosh(arg);
      }

      /// \copydoc UnaryJacobianOperation::getMathRep()
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

  /// UnaryJacobianOperation implementation for sqrt
  template<typename T_Real>
  struct OperationSqrt : public UnaryJacobianOperation<T_Real, OperationSqrt<T_Real>> {
    public:

      using Real = CODI_DD(T_Real, double);  ///< See UnaryJacobianOperation.

      /// \copydoc UnaryJacobianOperation::primal
      template<typename Arg>
      static CODI_INLINE Real primal(Arg const& arg) {
        return sqrt(arg);
      }

      /// \copydoc UnaryJacobianOperation::gradient
      template<typename Arg>
      static CODI_INLINE Real gradient(Arg const& arg, Real const& result) {
        if (Config::CheckExpressionArguments) {
          checkArgument(arg);
        }
        if (result != 0.0) {
          return 0.5 / result;
        } else {
          return (Real)0.0;
        }
      }

      /// \copydoc UnaryJacobianOperation::getMathRep()
      static CODI_INLINE std::string getMathRep() {
        return "sqrt";
      }

    private:
      template<typename Arg>
      CODI_INLINE static void checkArgument(Arg const& arg) {
        if (0.0 > RealTraits::getPassiveValue(arg)) {
          CODI_EXCEPTION("Sqrt of negative value or zero. (Value: %0.15e)", RealTraits::getPassiveValue(arg));
        }
      }

      template<typename Arg>
      CODI_INLINE static void checkArgument(std::complex<Arg> const& arg) {
        CODI_UNUSED(arg);
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

  /// UnaryJacobianOperation implementation for tan
  template<typename T_Real>
  struct OperationTan : public UnaryJacobianOperation<T_Real, OperationTan<T_Real>> {
    public:

      using Real = CODI_DD(T_Real, double);  ///< See UnaryJacobianOperation.

      /// \copydoc UnaryJacobianOperation::primal
      template<typename Arg>
      static CODI_INLINE Real primal(Arg const& arg) {
        return tan(arg);
      }

      /// \copydoc UnaryJacobianOperation::gradient
      template<typename Arg>
      static CODI_INLINE Real gradient(Arg const& arg, Real const& result) {
        CODI_UNUSED(result);
        if (Config::CheckExpressionArguments) {
          checkArgument(arg);
        }
        Real tmp = 1.0 / cos(arg);
        return tmp * tmp;
      }

      /// \copydoc UnaryJacobianOperation::getMathRep()
      static CODI_INLINE std::string getMathRep() {
        return "tan";
      }

    private:
      template<typename Arg>
      CODI_INLINE static void checkArgument(Arg const& arg) {
        if (0.0 == abs(cos(RealTraits::getPassiveValue(arg)))) {
          CODI_EXCEPTION("Tan evaluated at (0.5 + i) * PI. (Value: %0.15e)", RealTraits::getPassiveValue(arg));
        }
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

  /// UnaryJacobianOperation implementation for tanh
  template<typename T_Real>
  struct OperationTanh : public UnaryJacobianOperation<T_Real, OperationTanh<T_Real>> {
    public:

      using Real = CODI_DD(T_Real, double);  ///< See UnaryJacobianOperation.

      /// \copydoc UnaryJacobianOperation::primal
      template<typename Arg>
      static CODI_INLINE Real primal(Arg const& arg) {
        return tanh(arg);
      }

      /// \copydoc UnaryJacobianOperation::gradient
      template<typename Arg>
      static CODI_INLINE Real gradient(Arg const& arg, Real const& result) {
        CODI_UNUSED(arg);
        return 1.0 - result * result;
      }

      /// \copydoc UnaryJacobianOperation::getMathRep()
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

  /// UnaryJacobianOperation implementation for tgamma
  template<typename T_Real>
  struct OperationTgamma : public UnaryJacobianOperation<T_Real, OperationTgamma<T_Real>> {
    public:

      using Real = CODI_DD(T_Real, double);  ///< See UnaryJacobianOperation.

      /// \copydoc UnaryJacobianOperation::primal
      template<typename Arg>
      static CODI_INLINE Real primal(Arg const& arg) {
        return tgamma(arg);
      }

      /// \copydoc UnaryJacobianOperation::gradient
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

      /// \copydoc UnaryJacobianOperation::getMathRep()
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
