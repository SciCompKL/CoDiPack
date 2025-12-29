/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2026 Chair for Scientific Computing (SciComp), University of Kaiserslautern-Landau
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

#include <complex>

#include "../../config.h"
#include "../../misc/macros.hpp"
#include "../computeExpression.hpp"
#include "../expressionInterface.hpp"
#include "../real/allOperators.hpp"
#include "complexPredef.hpp"
#include "realToComplexCast.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  template<typename T_Tape, typename T_ParallelToolbox>
  struct ParallelActiveType;

  using std::abs;
  using std::arg;
  using std::conj;
  using std::imag;
  using std::norm;
  using std::polar;
  using std::proj;
  using std::real;

  /*******************************************************************************/
  /// @name Builtin binary operators
  /// @{

// Use the Logic from the real definition
#define OPERATION_LOGIC OperationAdd
#define FUNCTION operator+
#include "binaryMixedComplexAndRealOverloads.tpp"

#define OPERATION_LOGIC OperationSubstract
#define FUNCTION operator-
#include "binaryMixedComplexAndRealOverloads.tpp"

#define OPERATION_LOGIC OperationMultiply
#define FUNCTION operator*
#include "binaryMixedComplexAndRealOverloads.tpp"

#define OPERATION_LOGIC OperationDivide
#define FUNCTION operator/
#include "binaryMixedComplexAndRealOverloads.tpp"

  /// @}
  /*******************************************************************************/
  /// @name Standard math library binary operators
  /// @{

  /// BinaryJacobianOperation specialization for complex polar.
  template<typename T_ComplexReal>
  struct OperationComplexPolar : public BinaryOperation<T_ComplexReal, OperationComplexPolar<T_ComplexReal>> {
    public:

      using ComplexReal = CODI_DD(T_ComplexReal, std::complex<double>);  ///< See OperationComplexPolar.

      using Real = typename ComplexReal::value_type;  ///< Inner type of the complex value.

      /// \copydoc codi::BinaryJacobianOperation::primal()
      template<typename ArgA, typename ArgB>
      static CODI_INLINE ComplexReal primal(ArgA const& argA, ArgB const& argB) {
        return polar(argA, argB);
      }

      /// \copydoc codi::BinaryOperation::applyTangentArgA
      template<typename Tangent, typename ArgA, typename ArgB>
      static CODI_INLINE auto applyTangentArgA(Tangent const& tangent, ComplexReal const& result, ArgA const& argA,
                                               ArgB const& argB);
      // TODO: Implement

      /// \copydoc codi::BinaryOperation::applyAdjointArgA
      template<typename Adjoint, typename ArgA, typename ArgB>
      static CODI_INLINE Real applyAdjointArgA(Adjoint const& adjoint, ComplexReal const& result, ArgA const& argA,
                                               ArgB const& argB) {
        CODI_UNUSED(result, argA);
        return cos(argB) * std::real(adjoint) + sin(argB) * std::imag(adjoint);
      }

      /// \copydoc codi::BinaryOperation::applyTangentArgB
      template<typename Tangent, typename ArgA, typename ArgB>
      static CODI_INLINE auto applyTangentArgB(Tangent const& tangent, ComplexReal const& result, ArgA const& argA,
                                               ArgB const& argB);
      // TODO: Implement

      /// \copydoc codi::BinaryOperation::applyAdjointArgB
      template<typename Adjoint, typename ArgA, typename ArgB>
      static CODI_INLINE Real applyAdjointArgB(Adjoint const& adjoint, ComplexReal const& result, ArgA const& argA,
                                               ArgB const& argB) {
        CODI_UNUSED(argA, argB);

        return -std::imag(result) * std::real(adjoint) + std::real(result) * std::imag(adjoint);
      }

      /// \copydoc codi::BinaryJacobianOperation::getMathRep()
      static CODI_INLINE std::string getMathRep() {
        return "polar()";
      }
  };

#define FUNCTION polar
#define OPERATION_LOGIC OperationComplexPolar
#include "binaryRealToComplexOverloads.tpp"

  /// BinaryJacobianOperation specialization for complex pow.
  template<typename T_Real>
  struct OperationPow<std::complex<T_Real>>
      : public BinaryJacobianOperation<std::complex<T_Real>, OperationPow<std::complex<T_Real>>> {
    public:

      using ComplexReal = CODI_DD(std::complex<T_Real>, std::complex<double>);  ///< See BinaryJacobianOperation.

      /// \copydoc codi::ComputeOperation::primal()
      template<typename ArgA, typename ArgB>
      static CODI_INLINE ComplexReal primal(ArgA const& argA, ArgB const& argB) {
        return pow(argA, argB);
      }

      /// \copydoc codi::BinaryJacobianOperation::gradientA
      template<typename ArgA, typename ArgB>
      static CODI_INLINE ComplexReal gradientA(ArgA const& argA, ArgB const& argB, ComplexReal const& result) {
        CODI_UNUSED(result);

        return argB * pow(argA, argB - 1.0);
      }

      /// \copydoc codi::BinaryJacobianOperation::gradientB
      template<typename ArgA, typename ArgB>
      static CODI_INLINE ComplexReal gradientB(ArgA const& argA, ArgB const& argB, ComplexReal const& result) {
        CODI_UNUSED(argB);

        // Complex cast for argA, since the real log for negative numbers is not defined.
        return log(ComplexReal(argA)) * result;
      }

      /// \copydoc codi::BinaryJacobianOperation::getMathRep()
      static CODI_INLINE std::string getMathRep() {
        return "pow()";
      }
  };

#define OPERATION_LOGIC OperationPow
#define FUNCTION pow
#include "binaryMixedComplexAndRealOverloads.tpp"

  /// @}
  /*******************************************************************************/
  /// @name Builtin binary comparison operators
  /// @{

#define OPERATOR ==
#include "conditionalBinaryMixedComplexAndRealOverloads.tpp"

#define OPERATOR !=
#include "conditionalBinaryMixedComplexAndRealOverloads.tpp"

  /// @}
  /*******************************************************************************/
  /// @name Standard math library unary operators
  /// @{

  // Functions handled by the real definitions:
  // exp, log, log10, sqrt, sin, cos, tan, asin, acos, atan, sinh, cosh, tanh, asinh, acosh, atanh,
  // Unary operators handled by the real definitions:
  // operator+, operator-

  /// UnaryJacobianOperation implementation for complex abs
  template<typename T_Real>
  struct OperationComplexAbs : public UnaryJacobianOperation<T_Real, OperationComplexAbs<T_Real>> {
    public:

      using Real = CODI_DD(T_Real, double);  ///< See UnaryJacobianOperation.
      using Jacobian = std::complex<Real>;   ///< See UnaryJacobianOperation.

      /// \copydoc UnaryJacobianOperation::primal
      template<typename Arg>
      static CODI_INLINE Real primal(Arg const& arg) {
        return abs(arg);
      }

      /// \copydoc UnaryJacobianOperation::gradient
      template<typename Arg>
      static CODI_INLINE Jacobian gradient(Arg const& arg, Real const& result) {
        checkResult(result);
        if (result != 0.0) {
          return Jacobian(real(arg) / result, -imag(arg) / result);
        } else {
          return Jacobian(0.0);
        }
      }

      /// \copydoc codi::BinaryJacobianOperation::getMathRep()
      static CODI_INLINE std::string getMathRep() {
        return "abs";
      }

    private:
      static CODI_INLINE void checkResult(Real const& result) {
        if (Config::CheckExpressionArguments) {
          if (RealTraits::getPassiveValue(result) == 0.0) {
            CODI_EXCEPTION("Zero divisor for abs derivative.");
          }
        }
      }
  };

#define FUNCTION abs
#define OPERATION_LOGIC OperationComplexAbs
#include "unaryComplexToRealOverloads.tpp"

  /// UnaryJacobianOperation implementation for complex arg
  template<typename T_Real>
  struct OperationComplexArg : public UnaryJacobianOperation<T_Real, OperationComplexArg<T_Real>> {
    public:

      using Real = CODI_DD(T_Real, double);  ///< See UnaryJacobianOperation.
      using Jacobian = std::complex<Real>;   ///< See UnaryJacobianOperation.

      /// \copydoc UnaryJacobianOperation::primal
      template<typename Argument>
      static CODI_INLINE Real primal(Argument const& argument) {
        return arg(argument);
      }

      /// \copydoc UnaryJacobianOperation::gradient
      template<typename Argument>
      static CODI_INLINE Jacobian gradient(Argument const& argument, Real const& result) {
        CODI_UNUSED(result);

        Real divisor = real(argument) * real(argument) + imag(argument) * imag(argument);
        checkDivisor(divisor);
        divisor = 1.0 / divisor;

        return Jacobian(-imag(argument) * divisor, -real(argument) * divisor);
      }

      /// \copydoc codi::BinaryJacobianOperation::getMathRep()
      static CODI_INLINE std::string getMathRep() {
        return "arg";
      }

    private:
      static CODI_INLINE void checkDivisor(Real const& devisor) {
        if (Config::CheckExpressionArguments) {
          if (RealTraits::getPassiveValue(devisor) == 0.0) {
            CODI_EXCEPTION("Zero divisor for arg derivative.");
          }
        }
      }
  };

#define FUNCTION arg
#define OPERATION_LOGIC OperationComplexArg
#include "unaryComplexToRealOverloads.tpp"

  /// UnaryJacobianOperation implementation for complex conj.
  template<typename T_ComplexReal>
  struct OperationComplexConj : public UnaryOperation<T_ComplexReal, OperationComplexConj<T_ComplexReal>> {
    public:

      using ComplexReal = CODI_DD(T_ComplexReal, std::complex<double>);  ///< See UnaryJacobianOperation.

      /// \copydoc UnaryOperation::primal
      template<typename Arg>
      static CODI_INLINE ComplexReal primal(Arg const& arg) {
        return conj(arg);
      }

      /// \copydoc UnaryOperation::applyTangentArg
      template<typename Tangent, typename Arg>
      static CODI_INLINE auto applyTangentArg(Tangent const& tangent, ComplexReal const& result, Arg const& arg) {
        CODI_UNUSED(arg, result);

        return conj(tangent);
      }

      /// \copydoc UnaryOperation::applyAdjointArg
      template<typename Adjoint, typename Arg>
      static CODI_INLINE auto applyAdjointArg(Adjoint const& adjoint, ComplexReal const& result, Arg const& arg) {
        CODI_UNUSED(arg, result);

        return conj(adjoint);
      }

      /// \copydoc codi::BinaryJacobianOperation::getMathRep()
      static CODI_INLINE std::string getMathRep() {
        return "conj";
      }
  };

#define FUNCTION conj
#define OPERATION_LOGIC OperationComplexConj
#include "../real/unaryOverloads.tpp"

  /// UnaryJacobianOperation implementation for complex imag.
  template<typename T_Real>
  struct OperationComplexImag : public UnaryJacobianOperation<T_Real, OperationComplexImag<T_Real>> {
    public:

      using Real = CODI_DD(T_Real, double);  ///< See UnaryJacobianOperation.
      using Jacobian = std::complex<Real>;   ///< See UnaryJacobianOperation.

      /// \copydoc UnaryJacobianOperation::primal
      template<typename Arg>
      static CODI_INLINE Real primal(Arg const& arg) {
        return arg.imag();
      }

      /// \copydoc UnaryJacobianOperation::gradient
      template<typename Arg>
      static CODI_INLINE Jacobian gradient(Arg const& arg, Real const& result) {
        CODI_UNUSED(arg, result);

        return Jacobian(0.0, -1.0);
      }

      /// \copydoc codi::BinaryJacobianOperation::getMathRep()
      static CODI_INLINE std::string getMathRep() {
        return "imag";
      }
  };

#define FUNCTION imag
#define OPERATION_LOGIC OperationComplexImag
#include "unaryComplexToRealOverloads.tpp"

  /// UnaryJacobianOperation implementation for complex real.
  template<typename T_Real>
  struct OperationComplexNorm : public UnaryJacobianOperation<T_Real, OperationComplexNorm<T_Real>> {
    public:

      using Real = CODI_DD(T_Real, double);  ///< See UnaryJacobianOperation.
      using Jacobian = std::complex<Real>;   ///< See UnaryJacobianOperation.

      /// \copydoc UnaryJacobianOperation::primal
      template<typename Arg>
      static CODI_INLINE Real primal(Arg const& arg) {
        return norm(arg);
      }

      /// \copydoc UnaryJacobianOperation::gradient
      template<typename Arg>
      static CODI_INLINE Jacobian gradient(Arg const& arg, Real const& result) {
        CODI_UNUSED(result);

        return Jacobian(2.0 * real(arg), -2.0 * imag(arg));
      }

      /// \copydoc codi::BinaryJacobianOperation::getMathRep()
      static CODI_INLINE std::string getMathRep() {
        return "norm";
      }
  };

#define FUNCTION norm
#define OPERATION_LOGIC OperationComplexNorm
#include "unaryComplexToRealOverloads.tpp"

  /// UnaryJacobianOperation implementation for complex proj.
  template<typename T_ComplexReal>
  struct OperationComplexProj : public UnaryJacobianOperation<T_ComplexReal, OperationComplexProj<T_ComplexReal>> {
    public:

      using ComplexReal = CODI_DD(T_ComplexReal, std::complex<double>);  ///< See UnaryJacobianOperation.
      using Jacobian = double;                                           ///< See UnaryJacobianOperation.

      /// \copydoc UnaryJacobianOperation::primal
      template<typename Arg>
      static CODI_INLINE ComplexReal primal(Arg const& argument) {
        return proj(argument);
      }

      /// \copydoc UnaryJacobianOperation::gradient
      template<typename Arg>
      static CODI_INLINE Jacobian gradient(Arg const& argument, ComplexReal const& result) {
        CODI_UNUSED(argument, result);

        return 1.0;
      }

      /// \copydoc codi::BinaryJacobianOperation::getMathRep()
      static CODI_INLINE std::string getMathRep() {
        return "proj";
      }
  };

#define FUNCTION proj
#define OPERATION_LOGIC OperationComplexProj
#include "../real/unaryOverloads.tpp"

  /// UnaryJacobianOperation implementation for complex real.
  template<typename T_Real>
  struct OperationComplexReal : public UnaryJacobianOperation<T_Real, OperationComplexReal<T_Real>> {
    public:

      using Real = CODI_DD(T_Real, double);  ///< See UnaryJacobianOperation.
      using Jacobian = std::complex<Real>;   ///< See UnaryJacobianOperation.

      /// \copydoc UnaryJacobianOperation::primal
      template<typename Arg>
      static CODI_INLINE Real primal(Arg const& arg) {
        return arg.real();
      }

      /// \copydoc UnaryJacobianOperation::gradient
      template<typename Arg>
      static CODI_INLINE Jacobian gradient(Arg const& arg, Real const& result) {
        CODI_UNUSED(arg, result);

        return Jacobian(1.0, 0.0);
      }

      /// \copydoc codi::BinaryJacobianOperation::getMathRep()
      static CODI_INLINE std::string getMathRep() {
        return "real";
      }
  };

#define FUNCTION real
#define OPERATION_LOGIC OperationComplexReal
#include "unaryComplexToRealOverloads.tpp"

  /// @}
}

namespace std {
  /// isfinite implementation for complex numbers.
  template<typename Real>
  bool isfinite(std::complex<Real> arg) {
    return isfinite(arg.real()) && isfinite(arg.imag());
  }

  using codi::abs;
  using codi::imag;
  using codi::real;

}
