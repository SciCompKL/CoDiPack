/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2023 Chair for Scientific Computing (SciComp), University of Kaiserslautern-Landau
 * Homepage: http://www.scicomp.uni-kl.de
 * Contact:  Prof. Nicolas R. Gauger (codi@scicomp.uni-kl.de)
 *
 * Lead developers: Max Sagebaum, Johannes Blühdorn (SciComp, University of Kaiserslautern-Landau)
 *
 * This file is part of CoDiPack (http://www.scicomp.uni-kl.de/software/codi).
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
#include <utility>

#include "../../config.h"
#include "../../misc/exceptions.hpp"
#include "../../misc/macros.hpp"
#include "../../traits/realTraits.hpp"
#include "../binaryExpression.hpp"
#include "../constantExpression.hpp"
#include "../lhsExpressionInterface.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /*******************************************************************************/
  /// @name Builtin binary operators
  /// @{

  /// BinaryOperation implementation for operator +
  template<typename T_Real>
  struct OperationAdd : public BinaryOperation<T_Real> {
    public:

      using Real = CODI_DD(T_Real, double);  ///< See BinaryOperation.

      /// \copydoc codi::BinaryOperation::primal()
      template<typename ArgA, typename ArgB>
      static CODI_INLINE Real primal(ArgA const& argA, ArgB const& argB) {
        return argA + argB;
      }

      /// \copydoc codi::BinaryOperation::gradientA()
      template<typename ArgA, typename ArgB>
      static CODI_INLINE RealTraits::PassiveReal<Real> gradientA(ArgA const& argA, ArgB const& argB,
                                                                 Real const& result) {
        CODI_UNUSED(argA, argB, result);
        return 1.0;
      }

      /// \copydoc codi::BinaryOperation::gradientB()
      template<typename ArgA, typename ArgB>
      static CODI_INLINE RealTraits::PassiveReal<Real> gradientB(ArgA const& argA, ArgB const& argB,
                                                                 Real const& result) {
        CODI_UNUSED(argA, argB, result);
        return 1.0;
      }
  };
#define OPERATION_LOGIC OperationAdd
#define FUNCTION operator+
#include "binaryOverloads.tpp"

  /// BinaryOperation implementation for operator -
  template<typename T_Real>
  struct OperationSubstract : public BinaryOperation<T_Real> {
    public:

      using Real = CODI_DD(T_Real, double);  ///< See BinaryOperation.

      /// \copydoc codi::BinaryOperation::primal()
      template<typename ArgA, typename ArgB>
      static CODI_INLINE Real primal(ArgA const& argA, ArgB const& argB) {
        return argA - argB;
      }

      /// \copydoc codi::BinaryOperation::gradientA()
      template<typename ArgA, typename ArgB>
      static CODI_INLINE RealTraits::PassiveReal<Real> gradientA(ArgA const& argA, ArgB const& argB,
                                                                 Real const& result) {
        CODI_UNUSED(argA, argB, result);
        return 1.0;
      }

      /// \copydoc codi::BinaryOperation::gradientB()
      template<typename ArgA, typename ArgB>
      static CODI_INLINE RealTraits::PassiveReal<Real> gradientB(ArgA const& argA, ArgB const& argB,
                                                                 Real const& result) {
        CODI_UNUSED(argA, argB, result);
        return -1.0;
      }
  };
#define OPERATION_LOGIC OperationSubstract
#define FUNCTION operator-
#include "binaryOverloads.tpp"

  /// BinaryOperation implementation for operator *
  template<typename T_Real>
  struct OperationMultiply : public BinaryOperation<T_Real> {
    public:

      using Real = CODI_DD(T_Real, double);  ///< See BinaryOperation.

      /// \copydoc codi::BinaryOperation::primal()
      template<typename ArgA, typename ArgB>
      static CODI_INLINE Real primal(ArgA const& argA, ArgB const& argB) {
        return argA * argB;
      }

      /// \copydoc codi::BinaryOperation::gradientA()
      template<typename ArgA, typename ArgB>
      static CODI_INLINE ArgB const& gradientA(ArgA const& argA, ArgB const& argB, Real const& result) {
        CODI_UNUSED(argA, argB, result);
        return argB;
      }

      /// \copydoc codi::BinaryOperation::gradientB()
      template<typename ArgA, typename ArgB>
      static CODI_INLINE ArgA const& gradientB(ArgA const& argA, ArgB const& argB, Real const& result) {
        CODI_UNUSED(argA, argB, result);
        return argA;
      }
  };
#define OPERATION_LOGIC OperationMultiply
#define FUNCTION operator*
#include "binaryOverloads.tpp"

  /// BinaryOperation implementation for operator /
  template<typename T_Real>
  struct OperationDivide : public BinaryOperation<T_Real> {
    public:

      using Real = CODI_DD(T_Real, double);  ///< See BinaryOperation.

      /// \copydoc codi::BinaryOperation::primal()
      template<typename ArgA, typename ArgB>
      static CODI_INLINE Real primal(ArgA const& argA, ArgB const& argB) {
        return argA / argB;
      }

      /// \copydoc codi::BinaryOperation::gradientA()
      template<typename ArgA, typename ArgB>
      static CODI_INLINE Real gradientA(ArgA const& argA, ArgB const& argB, Real const& result) {
        CODI_UNUSED(argA, result);

        checkArguments(argB);
        return 1.0 / argB;
      }

      /// \copydoc codi::BinaryOperation::gradientB()
      template<typename ArgA, typename ArgB>
      static CODI_INLINE Real gradientB(ArgA const& argA, ArgB const& argB, Real const& result) {
        CODI_UNUSED(argA, result);

        checkArguments(argB);
        return -result / argB;
      }

    private:
      template<typename ArgB>
      static CODI_INLINE void checkArguments(ArgB const& argB) {
        if (Config::CheckExpressionArguments) {
          if (RealTraits::isTotalZero(RealTraits::getPassiveValue(argB))) {
            CODI_EXCEPTION("Division called with divisor of zero.");
          }
        }
      }
  };
#define OPERATION_LOGIC OperationDivide
#define FUNCTION operator/
#include "binaryOverloads.tpp"

  /// @}
  /*******************************************************************************/
  /// @name Standard math library binary operators
  /// @{

  using std::atan2;
  using std::copysign;
  using std::fmax;
  using std::fmin;
  using std::frexp;
  using std::hypot;
  using std::ldexp;
  using std::max;
  using std::min;
  using std::pow;
  using std::remainder;

  /// BinaryOperation implementation for atan2
  template<typename T_Real>
  struct OperationAtan2 : public BinaryOperation<T_Real> {
    public:

      using Real = CODI_DD(T_Real, double);  ///< See BinaryOperation.

      /// \copydoc codi::BinaryOperation::primal()
      template<typename ArgA, typename ArgB>
      static CODI_INLINE Real primal(ArgA const& argA, ArgB const& argB) {
        return atan2(argA, argB);
      }

      /// \copydoc codi::BinaryOperation::gradientA()
      template<typename ArgA, typename ArgB>
      static CODI_INLINE Real gradientA(ArgA const& argA, ArgB const& argB, Real const& result) {
        CODI_UNUSED(result);

        checkArguments(argA, argB);
        Real divisor = argA * argA + argB * argB;
        divisor = 1.0 / divisor;
        return argB * divisor;
      }

      /// \copydoc codi::BinaryOperation::gradientB()
      template<typename ArgA, typename ArgB>
      static CODI_INLINE Real gradientB(ArgA const& argA, ArgB const& argB, Real const& result) {
        CODI_UNUSED(result);

        checkArguments(argA, argB);
        Real divisor = argA * argA + argB * argB;
        divisor = 1.0 / divisor;
        return -argA * divisor;
      }

    private:
      template<typename ArgA, typename ArgB>
      static CODI_INLINE void checkArguments(ArgA& argA, ArgB& argB) {
        if (Config::CheckExpressionArguments) {
          if (0.0 == RealTraits::getPassiveValue(argA) && 0.0 == RealTraits::getPassiveValue(argB)) {
            CODI_EXCEPTION("atan2 called at point (0,0).");
          }
        }
      }
  };
#define OPERATION_LOGIC OperationAtan2
#define FUNCTION atan2
#include "binaryOverloads.tpp"

  /// BinaryOperation implementation for copysign
  template<typename T_Real>
  struct OperationCopysign : public BinaryOperation<T_Real> {
    public:

      using Real = CODI_DD(T_Real, double);  ///< See BinaryOperation.

      /// \copydoc codi::BinaryOperation::primal()
      template<typename ArgA, typename ArgB>
      static CODI_INLINE Real primal(ArgA const& argA, ArgB const& argB) {
        return copysign(argA, argB);
      }

      /// \copydoc codi::BinaryOperation::gradientA()
      template<typename ArgA, typename ArgB>
      static CODI_INLINE RealTraits::PassiveReal<Real> gradientA(ArgA const& argA, ArgB const& argB,
                                                                 Real const& result) {
        CODI_UNUSED(result);
        if (RealTraits::getPassiveValue(argA) < 0.0) {
          if (RealTraits::getPassiveValue(argB) < 0.0) {
            return RealTraits::PassiveReal<Real>(1.0);
          } else {
            return RealTraits::PassiveReal<Real>(-1.0);
          }
        } else if (RealTraits::getPassiveValue(argA) > 0.0) {
          if (RealTraits::getPassiveValue(argB) < 0.0) {
            return RealTraits::PassiveReal<Real>(-1.0);
          } else {
            return RealTraits::PassiveReal<Real>(1.0);
          }
        } else {
          return RealTraits::PassiveReal<Real>(0.0);
        }
      }

      /// \copydoc codi::BinaryOperation::gradientB()
      template<typename ArgA, typename ArgB>
      static CODI_INLINE RealTraits::PassiveReal<Real> gradientB(ArgA const& argA, ArgB const& argB,
                                                                 Real const& result) {
        CODI_UNUSED(argA, argB, result);

        return RealTraits::PassiveReal<Real>(0.0);
      }
  };

#define OPERATION_LOGIC OperationCopysign
#define FUNCTION copysign
#include "binaryOverloads.tpp"

#define OPERATION_LOGIC OperationCopysign
#define FUNCTION copysignf
#include "binaryOverloads.tpp"

  /// BinaryOperation implementation for frexp
  template<typename T_Real>
  struct OperationFrexp : public BinaryOperation<T_Real> {
    public:

      using Real = CODI_DD(T_Real, double);  ///< See BinaryOperation.

      /// \copydoc codi::BinaryOperation::primal()
      template<typename ArgA, typename ArgB>
      static CODI_INLINE Real primal(ArgA const& argA, ArgB const& argB) {
        return frexp(argA, argB);
      }

      /// \copydoc codi::BinaryOperation::gradientA()
      template<typename ArgA, typename ArgB>
      static CODI_INLINE Real gradientA(ArgA const& argA, ArgB const& argB, Real const& result) {
        CODI_UNUSED(argA, result);
        /* The result is always computed beforehand, therefore we can safely use the value of b. */
        return ldexp(1.0, -(*argB));
      }

      /// \copydoc codi::BinaryOperation::gradientB()
      template<typename ArgA, typename ArgB>
      static CODI_INLINE Real gradientB(ArgA const& argA, ArgB const& argB, Real const& result) {
        CODI_UNUSED(argA, argB, result);

        return Real();
      }
  };

#ifndef DOXYGEN_DISABLE
  template<typename T_StoreData>
  struct IntPointerConversion : public ConstantDataConversion<T_StoreData> {
    public:

      using StoreData = CODI_DD(T_StoreData, double);
      using ArgumentData = int*;

      static int* fromDataStore(StoreData const& v) {
        static int i = (int)v;
        return &i;
      }

      static StoreData toDataStore(int const* v) {
        return (StoreData)*v;
      }
  };
#endif
#define OPERATION_LOGIC OperationFrexp
#define FUNCTION frexp
#define SECOND_ARG_TYPE int*
#define SECOND_ARG_CONVERSION IntPointerConversion
#include "binaryFirstArgumentOverloads.tpp"

  /// BinaryOperation implementation for hypot
  template<typename T_Real>
  struct OperationHypot : public BinaryOperation<T_Real> {
    public:

      using Real = CODI_DD(T_Real, double);  ///< See BinaryOperation.

      /// \copydoc codi::BinaryOperation::primal()
      template<typename ArgA, typename ArgB>
      static CODI_INLINE Real primal(ArgA const& argA, ArgB const& argB) {
        return hypot(argA, argB);
      }

      /// \copydoc codi::BinaryOperation::gradientA()
      template<typename ArgA, typename ArgB>
      static CODI_INLINE Real gradientA(ArgA const& argA, ArgB const& argB, Real const& result) {
        CODI_UNUSED(argB);

        checkResult(result);
        if (result != 0.0) {
          return argA / result;
        } else {
          return Real();
        }
      }

      /// \copydoc codi::BinaryOperation::gradientB()
      template<typename ArgA, typename ArgB>
      static CODI_INLINE Real gradientB(ArgA const& argA, ArgB const& argB, Real const& result) {
        CODI_UNUSED(argA);

        checkResult(result);
        if (result != 0.0) {
          return argB / result;
        } else {
          return Real();
        }
      }

    private:
      static CODI_INLINE void checkResult(Real const& result) {
        if (Config::CheckExpressionArguments) {
          if (RealTraits::getPassiveValue(result) == 0.0) {
            CODI_EXCEPTION("Zero divisor for hypot derivative.");
          }
        }
      }
  };

#define OPERATION_LOGIC OperationHypot
#define FUNCTION hypot
#include "binaryOverloads.tpp"

#define OPERATION_LOGIC OperationHypot
#define FUNCTION hypotl
#include "binaryOverloads.tpp"

#define OPERATION_LOGIC OperationHypot
#define FUNCTION hypotf
#include "binaryOverloads.tpp"

  /// BinaryOperation implementation for ldexp
  template<typename T_Real>
  struct OperationLdexp : public BinaryOperation<T_Real> {
    public:

      using Real = CODI_DD(T_Real, double);  ///< See BinaryOperation.

      /// \copydoc codi::BinaryOperation::primal()
      template<typename ArgA, typename ArgB>
      static CODI_INLINE Real primal(ArgA const& argA, ArgB const& argB) {
        return ldexp(argA, argB);
      }

      /// \copydoc codi::BinaryOperation::gradientA()
      template<typename ArgA, typename ArgB>
      static CODI_INLINE Real gradientA(ArgA const& argA, ArgB const& argB, Real const& result) {
        CODI_UNUSED(argA, result);

        return ldexp(1.0, argB);
      }

      /// \copydoc codi::BinaryOperation::gradientB()
      template<typename ArgA, typename ArgB>
      static CODI_INLINE Real gradientB(ArgA const& argA, ArgB const& argB, Real const& result) {
        CODI_UNUSED(argA, argB, result);

        return Real();
      }
  };
#define OPERATION_LOGIC OperationLdexp
#define FUNCTION ldexp
#define SECOND_ARG_TYPE int
#define SECOND_ARG_CONVERSION ConstantDataConversion
#include "binaryFirstArgumentOverloads.tpp"

  /// BinaryOperation implementation for max
  template<typename T_Real>
  struct OperationMax : public BinaryOperation<T_Real> {
    public:

      using Real = CODI_DD(T_Real, double);  ///< See BinaryOperation.

      /// \copydoc codi::BinaryOperation::primal()
      template<typename ArgA, typename ArgB>
      static CODI_INLINE Real primal(ArgA const& argA, ArgB const& argB) {
        return max(argA, argB);
      }

      /// \copydoc codi::BinaryOperation::gradientA()
      template<typename ArgA, typename ArgB>
      static CODI_INLINE RealTraits::PassiveReal<Real> gradientA(ArgA const& argA, ArgB const& argB,
                                                                 Real const& result) {
        CODI_UNUSED(result);

        if (RealTraits::getPassiveValue(argA) > RealTraits::getPassiveValue(argB)) {
          return 1.0;
        } else {
          return 0.0;
        }
      }

      /// \copydoc codi::BinaryOperation::gradientB()
      template<typename ArgA, typename ArgB>
      static CODI_INLINE RealTraits::PassiveReal<Real> gradientB(ArgA const& argA, ArgB const& argB,
                                                                 Real const& result) {
        CODI_UNUSED(result);

        if (RealTraits::getPassiveValue(argA) > RealTraits::getPassiveValue(argB)) {
          return 0.0;
        } else {
          return 1.0;
        }
      }
  };

#define OPERATION_LOGIC OperationMax
#define FUNCTION max
#include "binaryOverloads.tpp"

#define OPERATION_LOGIC OperationMax
#define FUNCTION fmax
#include "binaryOverloads.tpp"

  /// BinaryOperation implementation for min
  template<typename T_Real>
  struct OperationMin : public BinaryOperation<T_Real> {
    public:

      using Real = CODI_DD(T_Real, double);  ///< See BinaryOperation.

      /// \copydoc codi::BinaryOperation::primal()
      template<typename ArgA, typename ArgB>
      static CODI_INLINE Real primal(ArgA const& argA, ArgB const& argB) {
        return min(argA, argB);
      }

      /// \copydoc codi::BinaryOperation::gradientA()
      template<typename ArgA, typename ArgB>
      static CODI_INLINE RealTraits::PassiveReal<Real> gradientA(ArgA const& argA, ArgB const& argB,
                                                                 Real const& result) {
        CODI_UNUSED(result);

        if (RealTraits::getPassiveValue(argA) < RealTraits::getPassiveValue(argB)) {
          return 1.0;
        } else {
          return 0.0;
        }
      }

      /// \copydoc codi::BinaryOperation::gradientB()
      template<typename ArgA, typename ArgB>
      static CODI_INLINE RealTraits::PassiveReal<Real> gradientB(ArgA const& argA, ArgB const& argB,
                                                                 Real const& result) {
        CODI_UNUSED(result);

        if (RealTraits::getPassiveValue(argA) < RealTraits::getPassiveValue(argB)) {
          return 0.0;
        } else {
          return 1.0;
        }
      }
  };
#define OPERATION_LOGIC OperationMin
#define FUNCTION min
#include "binaryOverloads.tpp"

#define OPERATION_LOGIC OperationMin
#define FUNCTION fmin
#include "binaryOverloads.tpp"

  /// BinaryOperation implementation for pow
  template<typename T_Real>
  struct OperationPow : public BinaryOperation<T_Real> {
    public:

      using Real = CODI_DD(T_Real, double);  ///< See BinaryOperation.

      /// \copydoc codi::BinaryOperation::primal()
      template<typename ArgA, typename ArgB>
      static CODI_INLINE Real primal(ArgA const& argA, ArgB const& argB) {
        return pow(argA, argB);
      }

      /// \copydoc codi::BinaryOperation::gradientA()
      template<typename ArgA, typename ArgB>
      static CODI_INLINE Real gradientA(ArgA const& argA, ArgB const& argB, Real const& result) {
        CODI_UNUSED(result);

        checkArguments(argA);
        if (RealTraits::getPassiveValue(argA) <= 0.0 && 1 <= RealTraits::MaxDerivativeOrder<ArgB>()) {
          // Special case for higher order derivatives. Derivative will be wrong since the argB part is not evaluated.
          return RealTraits::getPassiveValue(argB) * pow(argA, argB - 1.0);
        } else {
          return argB * pow(argA, argB - 1.0);
        }
      }

      /// \copydoc codi::BinaryOperation::gradientB()
      template<typename ArgA, typename ArgB>
      static CODI_INLINE Real gradientB(ArgA const& argA, ArgB const& argB, Real const& result) {
        CODI_UNUSED(argB);

        checkArguments(argA);
        if (RealTraits::getPassiveValue(argA) > 0.0) {
          return log(argA) * result;
        } else {
          return Real();
        }
      }

    private:
      template<typename ArgA>
      static CODI_INLINE void checkArguments(ArgA& argA) {
        if (Config::CheckExpressionArguments) {
          if (RealTraits::getPassiveValue(argA) < 0.0) {
            CODI_EXCEPTION("Negative base for active exponent in pow function. (Value: %0.15e)",
                           RealTraits::getPassiveValue(argA));
          }
        }
      }
  };
#define OPERATION_LOGIC OperationPow
#define FUNCTION pow
#include "binaryOverloads.tpp"

  /// BinaryOperation implementation for remainder
  ///
  /// Derivative implementation based on IEC 60559: remainder = numer - rquot * denom
  template<typename T_Real>
  struct OperationRemainder : public BinaryOperation<T_Real> {
    public:

      using Real = CODI_DD(T_Real, double);  ///< See BinaryOperation.

      /// \copydoc codi::BinaryOperation::primal()
      template<typename ArgA, typename ArgB>
      static CODI_INLINE Real primal(ArgA const& argA, ArgB const& argB) {
        return remainder(argA, argB);
      }

      /// \copydoc codi::BinaryOperation::gradientA()
      template<typename ArgA, typename ArgB>
      static CODI_INLINE Real gradientA(ArgA const& argA, ArgB const& argB, Real const& result) {
        CODI_UNUSED(argA, argB, result);

        return (Real)1.0;
      }

      /// \copydoc codi::BinaryOperation::gradientB()
      template<typename ArgA, typename ArgB>
      static CODI_INLINE Real gradientB(ArgA const& argA, ArgB const& argB, Real const& result) {
        CODI_UNUSED(result);

        checkArguments(argB);

        using std::round;
        return -round(argA / argB);
      }

    private:
      template<typename ArgB>
      static CODI_INLINE void checkArguments(ArgB const& argB) {
        if (Config::CheckExpressionArguments) {
          if (0.0 == RealTraits::getPassiveValue(argB)) {
            CODI_EXCEPTION("Remainder called with divisor of zero.");
          }
        }
      }
  };
#define OPERATION_LOGIC OperationRemainder
#define FUNCTION remainder
#include "binaryOverloads.tpp"

  /// @}
  /// /*******************************************************************************/
  /// @name Additional standard library binary operators
  /// @{

  using std::swap;

  /// Optimized swap for lhs expressions that does not call the copy constructor.
  template<typename Real, typename Gradient, typename Tape, typename Impl>
  void swap(LhsExpressionInterface<Real, Gradient, Tape, Impl>& lhs,
            LhsExpressionInterface<Real, Gradient, Tape, Impl>& rhs) {
    swap(lhs.value(), rhs.value());
    swap(lhs.getIdentifier(), rhs.getIdentifier());
  }

  /// @}
}

namespace std {
  using codi::atan2;
  using codi::copysign;
  using codi::copysignf;
  using codi::fmax;
  using codi::fmin;
  using codi::frexp;
  using codi::hypot;
  using codi::hypotf;
  using codi::hypotl;
  using codi::ldexp;
  using codi::max;
  using codi::min;
  using codi::pow;
  using codi::remainder;
  using codi::swap;
}
