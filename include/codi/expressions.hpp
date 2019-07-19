/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2019 Chair for Scientific Computing (SciComp), TU Kaiserslautern
 * Homepage: http://www.scicomp.uni-kl.de
 * Contact:  Prof. Nicolas R. Gauger (codi@scicomp.uni-kl.de)
 *
 * Lead developers: Max Sagebaum, Tim Albring (SciComp, TU Kaiserslautern)
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
 * Authors: Max Sagebaum, Tim Albring, (SciComp, TU Kaiserslautern)
 *          Prof. Robin Hogan, (Univ. of Reading).
 *
 * Originally based on Adept 1.0 (http://www.met.rdg.ac.uk/clouds/adept/)
 * released under GPL 3.0 (Copyright (C) 2012-2013 Robin Hogan and the University of Reading).
 */

#pragma once

#include <cmath>
#include <algorithm>
#include <iostream>

#include "configure.h"
#include "exceptions.hpp"
#include "macros.h"
#include "typeTraits.hpp"

/**
 * @brief Global namespace for CoDiPack - Code Differentiation Package
 */
namespace codi {

  template<typename A> struct ExpressionTraits;

  /**
   * The Expression type from which all other types of expression
   * derive. Each member function simply calls the specialized version
   * of the function according to the expression's true type, which is
   * given by its template argument.
   *
   * @tparam Real  The data type of the primal values and the gradient values.
   * @tparam    A  The implementing class of the expression.
   */
  template<typename Real, class A>
  struct Expression {

    /**
     * @brief If true, implementations of the expression are stored as references otherwise by values.
     *
     * This values is used by the macro CODI_CREATE_STORE_TYPE.
     */
    static const bool storeAsReference;

    /**
     * @brief The passive value is used where the expressions are combined with normal double values.
     */
    typedef typename TypeTraits<Real>::PassiveReal PassiveReal;

    /**
     * Cast the expression to its true type, given by the template
     * argument.
     *
     * @return The instance of the implementing class.
     */
    CODI_INLINE const A& cast() const {
      return static_cast<const A&>(*this);
    }

    /**
     * @brief Calculate the gradient of the expression.
     *
     * Calculate the gradient of the mathematical operation that this
     * expression represents and pass the result to its argument.
     * For functions f(a), pass df/da to the argument.
     *
     * @param[in,out] data A helper value which the tape can define and use for the evaluation.
     *
     * @tparam Data The type for the tape data.
     */
    template<typename Data>
    CODI_INLINE void calcGradient(Data& data) const {
      cast().calcGradient(data);
    }

    /**
     * @brief Calculate the gradient of the expression.
     *
     * Calculate the gradient of the mathematical operation that this
     * expression represents and pass the result to its argument.
     * For functions f(a), pass multiplier * df/da to the argument.
     *
     * @param[in,out]    data A helper value which the tape can define and use for the evaluation.
     * @param[in]  multiplier The Jacobi from the expression where this expression was used as an argument.
     *
     * @tparam Data The type for the tape data.
     */
    template<typename Data>
    CODI_INLINE void calcGradient(Data& data, const Real& multiplier) const {
      cast().calcGradient(data, multiplier);
    }

    /**
     * @brief Return the numerical value of the expression.
     *
     * @return The value of the expression.
     */
    CODI_INLINE const Real getValue() const {
      return cast().getValue();
    }

    /**
     * @brief constantValueActions are called for every constant real in the expression.
     *
     * @param[in,out] tape  The tape that calls the action.
     * @param[in,out] data  The data that can be used by the action.
     * @param[in]     func  The function that is called for every constant item.
     *
     * @tparam CallTape  The type of the tape that calls the action.
     * @tparam     Data  The type of the data for the action.
     * @tparam     Func  The type of the function that is called.
     */
    template<typename CallTape, typename Data, typename Func>
    CODI_INLINE void constantValueAction(CallTape& tape, Data data, Func func) const {
      cast().constantValueAction(tape, data, func);
    }

    /**
     * @brief The action is called for every active real in the expression.
     *
     * @param[in,out] data  The data that can be used by the action.
     * @param[in]     func  The function that is called for every active real in the expression.
     *
     * @tparam     Data  The type of the data for the action.
     * @tparam     Func  The type of the function that is called.
     */
    template<typename Data, typename Func>
    CODI_INLINE void valueAction(Data data, Func func) const {
      cast().valueAction(data, func);
    }

#if CODI_EnableImplicitConversion
    /**
     * @brief Get the primal value of this instance via implicit cast.
     * @return The primal value.
     */
    CODI_INLINE operator const Real() const {
      Warning::implicitCast<CODI_DisableImplicitConversionWarning>();

      return getValue();
    }
#endif

  private:
    /**
     * Intentionally inaccessible to prevent an expression appearing
     * on the left-hand-side of a statement
     */
    Expression& operator=(const Expression&) = delete;
  };

  /*
   * Now define particular types of expression, using static
   * polymorphism via the Curiously Recurring Template Pattern
   */

  /**
   * @brief The macro creates a helper function that calls an operator as a function.
   *
   * The generated function has the format
   *
   * CODI_INLINE auto primal_NAME(const A& a, const B& b);
   *
   * @param   NAME  The name for the generated function.
   * @param   OP    The sign of the operator that the function calls.
   */
  #define CODI_OPERATOR_HELPER(NAME, OP) \
    /** @brief Helper function to call operators as a function. @param[in] a The first argument of the operation. @param[in] b The second argument of the operation. @return The value of a OP b @tparam A The expression for the first argument of the function @tparam B The expression for the second argument of the function*/ \
    template<typename A, typename B> \
    CODI_INLINE auto primal_ ## NAME(const A& a, const B& b) -> decltype(a OP b) { \
      return a OP b; \
    }

  /*
   * Now all binary functions are implemented.
   *
   * We use the naming scheme dervBB[M]_Name.
   *
   * BB tells which variable is active. Thus we have the combinations
   * 11 -> both are active
   * 10 -> first argument is active
   * 01 -> second argument is active
   *
   * There are also the functions gradientA_Name and gradientB_Name which calculate the jacobie with respect to the first and
   * second argument respectively.
   *
   * There is no implementation for 00 because no variable is active and thus the derivative would be zero.
   *
   * If the M is present the method is implemented with the multiplier as an argument.
   */

  /*
   * Implementation for f(a,b) = a + b
   * df/da = 1 and likewise for df/db so simply
   * call a and b's versions of calcGradient
   */
  template<typename Real, typename A, typename B> CODI_INLINE const typename TypeTraits<Real>::PassiveReal gradientA_Add(const A& a, const B& b, const Real& result) {
    CODI_UNUSED(a);
    CODI_UNUSED(b);
    CODI_UNUSED(result);
    return 1.0;
  }
  template<typename Real, typename A, typename B> CODI_INLINE const typename TypeTraits<Real>::PassiveReal gradientB_Add(const A& a, const B& b, const Real& result) {
    CODI_UNUSED(a);
    CODI_UNUSED(b);
    CODI_UNUSED(result);
    return 1.0;
  }
  template<typename Data, typename Real, typename A, typename B> CODI_INLINE void derv11_Add(Data& data, const A& a, const B& b, const Real& result) {
    CODI_UNUSED(result);
    a.calcGradient(data);
    b.calcGradient(data);
  }
  template<typename Data, typename Real, typename A, typename B> CODI_INLINE void derv11M_Add(Data& data, const A& a, const B& b, const Real& result, const Real& multiplier) {
    CODI_UNUSED(result);
    a.calcGradient(data, multiplier);
    b.calcGradient(data, multiplier);
  }
  template<typename Data, typename Real, typename A> CODI_INLINE void derv10_Add(Data& data, const A& a, const typename TypeTraits<Real>::PassiveReal& b, const Real& result) {
    CODI_UNUSED(b);
    CODI_UNUSED(result);
    a.calcGradient(data);
  }
  template<typename Data, typename Real, typename A> CODI_INLINE void derv10M_Add(Data& data, const A& a, const typename TypeTraits<Real>::PassiveReal& b, const Real& result, const Real& multiplier) {
    CODI_UNUSED(b);
    CODI_UNUSED(result);
    a.calcGradient(data, multiplier);
  }
  template<typename Data, typename Real, typename B> CODI_INLINE void derv01_Add(Data& data, const typename TypeTraits<Real>::PassiveReal& a, const B& b, const Real& result) {
    CODI_UNUSED(a);
    CODI_UNUSED(result);
    b.calcGradient(data);
  }
  template<typename Data, typename Real, typename B> CODI_INLINE void derv01M_Add(Data& data, const typename TypeTraits<Real>::PassiveReal& a, const B& b, const Real& result, const Real& multiplier) {
    CODI_UNUSED(a);
    CODI_UNUSED(result);
    b.calcGradient(data, multiplier);
  }

  CODI_OPERATOR_HELPER(Add, +)
  #define NAME Add
  #define FUNCTION operator +
  #define PRIMAL_FUNCTION primal_Add
  #include "binaryExpression.tpp"

  /*
   * Implementation for f(a,b) = a - b
   * df/da = 1 so simply call a
   */
  template<typename Real, typename A, typename B> CODI_INLINE const typename TypeTraits<Real>::PassiveReal gradientA_Subtract(const A& a, const B& b, const Real& result) {
    CODI_UNUSED(a);
    CODI_UNUSED(b);
    CODI_UNUSED(result);
    return 1.0;
  }
  template<typename Real, typename A, typename B> CODI_INLINE const typename TypeTraits<Real>::PassiveReal gradientB_Subtract(const A& a, const B& b, const Real& result) {
    CODI_UNUSED(a);
    CODI_UNUSED(b);
    CODI_UNUSED(result);
    return -1.0;
  }
  template<typename Data, typename Real, typename A, typename B> CODI_INLINE void derv11_Subtract(Data& data, const A& a, const B& b, const Real& result) {
    CODI_UNUSED(result);
    a.calcGradient(data);
    b.calcGradient(data, -1.0);
  }
  template<typename Data, typename Real, typename A, typename B> CODI_INLINE void derv11M_Subtract(Data& data, const A& a, const B& b, const Real& result, const Real& multiplier) {
    CODI_UNUSED(result);
    a.calcGradient(data, multiplier);
    b.calcGradient(data, -multiplier);
  }
  template<typename Data, typename Real, typename A> CODI_INLINE void derv10_Subtract(Data& data, const A& a, const typename TypeTraits<Real>::PassiveReal& b, const Real& result) {
    CODI_UNUSED(b);
    CODI_UNUSED(result);
    a.calcGradient(data);
  }
  template<typename Data, typename Real, typename A> CODI_INLINE void derv10M_Subtract(Data& data, const A& a, const typename TypeTraits<Real>::PassiveReal& b, const Real& result, const Real& multiplier) {
    CODI_UNUSED(b);
    CODI_UNUSED(result);
    a.calcGradient(data, multiplier);
  }
  template<typename Data, typename Real, typename B> CODI_INLINE void derv01_Subtract(Data& data, const typename TypeTraits<Real>::PassiveReal& a, const B& b, const Real& result) {
    CODI_UNUSED(a);
    CODI_UNUSED(result);
    b.calcGradient(data, -1.0);
  }
  template<typename Data, typename Real, typename B> CODI_INLINE void derv01M_Subtract(Data& data, const typename TypeTraits<Real>::PassiveReal& a, const B& b, const Real& result, const Real& multiplier) {
    CODI_UNUSED(a);
    CODI_UNUSED(result);
    b.calcGradient(data, -multiplier);
  }
  CODI_OPERATOR_HELPER(Subtract, -)
  #define NAME Subtract
  #define FUNCTION operator -
  #define PRIMAL_FUNCTION primal_Subtract
  #include "binaryExpression.tpp"

  /*
   * Implementation for f(a,b) = a * b
   */
  template<typename Real, typename A, typename B> CODI_INLINE const B& gradientA_Multiply(const A& a, const B& b, const Real& result) {
    CODI_UNUSED(a);
    CODI_UNUSED(result);
    return b;
  }
  template<typename Real, typename A, typename B> CODI_INLINE const A& gradientB_Multiply(const A& a, const B& b, const Real& result) {
    CODI_UNUSED(b);
    CODI_UNUSED(result);
    return a;
  }
  template<typename Data, typename Real, typename A, typename B> CODI_INLINE void derv11_Multiply(Data& data, const A& a, const B& b, const Real& result) {
    CODI_UNUSED(result);
    a.calcGradient(data, b.getValue());
    b.calcGradient(data, a.getValue());
  }
  template<typename Data, typename Real, typename A, typename B> CODI_INLINE void derv11M_Multiply(Data& data, const A& a, const B& b, const Real& result, const Real& multiplier) {
    CODI_UNUSED(result);
    a.calcGradient(data, b.getValue() * multiplier);
    b.calcGradient(data, a.getValue() * multiplier);
  }
  template<typename Data, typename Real, typename A> CODI_INLINE void derv10_Multiply(Data& data, const A& a, const typename TypeTraits<Real>::PassiveReal& b, const Real& result) {
    CODI_UNUSED(result);
    a.calcGradient(data, b);
  }
  template<typename Data, typename Real, typename A> CODI_INLINE void derv10M_Multiply(Data& data, const A& a, const typename TypeTraits<Real>::PassiveReal& b, const Real& result, const Real& multiplier) {
    CODI_UNUSED(result);
    a.calcGradient(data, b * multiplier);
  }
  template<typename Data, typename Real, typename B> CODI_INLINE void derv01_Multiply(Data& data, const typename TypeTraits<Real>::PassiveReal& a, const B& b, const Real& result) {
    CODI_UNUSED(result);
    b.calcGradient(data, a);
  }
  template<typename Data, typename Real, typename B> CODI_INLINE void derv01M_Multiply(Data& data, const typename TypeTraits<Real>::PassiveReal& a, const B& b, const Real& result, const Real& multiplier) {
    CODI_UNUSED(result);
    b.calcGradient(data, a * multiplier);
  }
  CODI_OPERATOR_HELPER(Multiply, *)
  #define NAME Multiply
  #define FUNCTION operator *
  #define PRIMAL_FUNCTION primal_Multiply
  #include "binaryExpression.tpp"

  /*
   * Implementation for f(a,b) = a / b
   */
  /**
   * @brief Helper function which checks if the divisor is zero.
   *
   * Depending of the global option CheckExpressionArguments the function will call
   * CODI_EXCEPTION if b is equal to zero.
   *
   * @tparam Real The real type used in the active type.
   */
  template<typename Real> CODI_INLINE void checkArgumentsDivide(const Real& b) {
    if(CheckExpressionArguments) {
      if( 0.0 == TypeTraits<Real>::getBaseValue(b)) {
        CODI_EXCEPTION("Division called with divisor of zero.");
      }
    }
  }
  template<typename Real, typename A, typename B> CODI_INLINE const Real gradientA_Divide(const A& a, const B& b, const Real& result) {
    checkArgumentsDivide(b);
    CODI_UNUSED(a);
    CODI_UNUSED(result);
    return 1.0 / b;
  }
  template<typename Real, typename A, typename B> CODI_INLINE const Real gradientB_Divide(const A& a, const B& b, const Real& result) {
    checkArgumentsDivide(b);
    CODI_UNUSED(a);
    return -result / b;
  }
  template<typename Data, typename Real, typename A, typename B> CODI_INLINE void derv11_Divide(Data& data, const A& a, const B& b, const Real& result) {
    checkArgumentsDivide(b.getValue());
    Real one_over_b = 1.0 / b.getValue();
    a.calcGradient(data, one_over_b);
    b.calcGradient(data, -result * one_over_b);
  }
  template<typename Data, typename Real, typename A, typename B> CODI_INLINE void derv11M_Divide(Data& data, const A& a, const B& b, const Real& result, const Real& multiplier) {
    checkArgumentsDivide(b.getValue());
    Real one_over_b = multiplier / b.getValue();
    a.calcGradient(data, one_over_b);
    b.calcGradient(data, -result * one_over_b);
  }
  template<typename Data, typename Real, typename A> CODI_INLINE void derv10_Divide(Data& data, const A& a, const typename TypeTraits<Real>::PassiveReal& b, const Real& result) {
    CODI_UNUSED(result);
    checkArgumentsDivide(b);
    Real one_over_b = 1.0 / b;
    a.calcGradient(data, one_over_b);
  }
  template<typename Data, typename Real, typename A> CODI_INLINE void derv10M_Divide(Data& data, const A& a, const typename TypeTraits<Real>::PassiveReal& b, const Real& result, const Real& multiplier) {
    CODI_UNUSED(result);
    checkArgumentsDivide(b);
    Real one_over_b = multiplier / b;
    a.calcGradient(data, one_over_b);
  }
  template<typename Data, typename Real, typename B> CODI_INLINE void derv01_Divide(Data& data, const typename TypeTraits<Real>::PassiveReal& a, const B& b, const Real& result) {
    CODI_UNUSED(a);
    checkArgumentsDivide(b.getValue());
    Real one_over_b = 1.0 / b.getValue();
    b.calcGradient(data, -result * one_over_b);
  }
  template<typename Data, typename Real, typename B> CODI_INLINE void derv01M_Divide(Data& data, const typename TypeTraits<Real>::PassiveReal& a, const B& b, const Real& result, const Real& multiplier) {
    CODI_UNUSED(a);
    checkArgumentsDivide(b.getValue());
    Real one_over_b = multiplier / b.getValue();
    b.calcGradient(data, -result * one_over_b);
  }
  CODI_OPERATOR_HELPER(Divide, /)
  #define NAME Divide
  #define FUNCTION operator /
  #define PRIMAL_FUNCTION primal_Divide
  #include "binaryExpression.tpp"

  /*
   * Implementation for f(a,b) = atan2(a,b)
   */
  /**
   * @brief Helper function which checks if both arguments are zero.
   *
   * Depending of the global option CheckExpressionArguments the function will call
   * CODI_EXCEPTION if a and b are equal to zero.
   *
   * @tparam A Type of the first argument.
   * @tparam B Type of the second argument.
   */
  template<typename A, typename B> CODI_INLINE void checkArgumentsAtan2(const A& a, const B& b) {
    if(CheckExpressionArguments) {
      if( 0.0 == TypeTraits<A>::getBaseValue(a) &&
          0.0 == TypeTraits<B>::getBaseValue(b)) {
        CODI_EXCEPTION("atan2 called at point (0,0).");
      }
    }
  }
  template<typename Real, typename A, typename B> CODI_INLINE const Real gradientA_Atan2(const A& a, const B& b, const Real& result) {
    CODI_UNUSED(result);
    checkArgumentsAtan2(a, b);
    Real divisor = a * a + b * b;
    divisor = 1.0 / divisor;
    return b * divisor;

  }
  template<typename Real, typename A, typename B> CODI_INLINE const Real gradientB_Atan2(const A& a, const B& b, const Real& result) {
    CODI_UNUSED(result);
    checkArgumentsAtan2(a, b);
    Real divisor = a * a + b * b;
    divisor = 1.0 / divisor;
    return -a * divisor;
  }
  template<typename Data, typename Real, typename A, typename B> CODI_INLINE void derv11_Atan2(Data& data, const A& a, const B& b, const Real& result) {
    CODI_UNUSED(result);
    checkArgumentsAtan2(a.getValue(), b.getValue());
    Real divisor = a.getValue() * a.getValue() + b.getValue() * b.getValue();
    divisor = 1.0 / divisor;
    a.calcGradient(data, b.getValue() * divisor);
    b.calcGradient(data, -a.getValue() * divisor);
  }
  template<typename Data, typename Real, typename A, typename B> CODI_INLINE void derv11M_Atan2(Data& data, const A& a, const B& b, const Real& result, const Real& multiplier) {
    CODI_UNUSED(result);
    checkArgumentsAtan2(a.getValue(), b.getValue());
    Real divisor = a.getValue() * a.getValue() + b.getValue() * b.getValue();
    divisor = 1.0 / divisor;
    a.calcGradient(data, multiplier * b.getValue() * divisor);
    b.calcGradient(data, multiplier * -a.getValue() * divisor);
  }
  template<typename Data, typename Real, typename A> CODI_INLINE void derv10_Atan2(Data& data, const A& a, const typename TypeTraits<Real>::PassiveReal& b, const Real& result) {
    CODI_UNUSED(result);
    checkArgumentsAtan2(a.getValue(), b);
    Real divisor = a.getValue() * a.getValue() + b * b;
    divisor = 1.0 / divisor;
    a.calcGradient(data, b * divisor);
  }
  template<typename Data, typename Real, typename A> CODI_INLINE void derv10M_Atan2(Data& data, const A& a, const typename TypeTraits<Real>::PassiveReal& b, const Real& result, const Real& multiplier) {
    CODI_UNUSED(result);
    checkArgumentsAtan2(a.getValue(), b);
    Real divisor = a.getValue() * a.getValue() + b * b;
    divisor = 1.0 / divisor;
    a.calcGradient(data, multiplier * b * divisor);
  }
  template<typename Data, typename Real, typename B> CODI_INLINE void derv01_Atan2(Data& data, const typename TypeTraits<Real>::PassiveReal& a, const B& b, const Real& result) {
    CODI_UNUSED(result);
    checkArgumentsAtan2(a, b.getValue());
    Real divisor = a * a + b.getValue() * b.getValue();
    divisor = 1.0 / divisor;
    b.calcGradient(data, -a * divisor);
  }
  template<typename Data, typename Real, typename B> CODI_INLINE void derv01M_Atan2(Data& data, const typename TypeTraits<Real>::PassiveReal& a, const B& b, const Real& result, const Real& multiplier) {
    CODI_UNUSED(result);
    checkArgumentsAtan2(a, b.getValue());
    Real divisor = a * a + b.getValue() * b.getValue();
    divisor = 1.0 / divisor;
    b.calcGradient(data, multiplier * -a * divisor);
  }
  using std::atan2;
  #define NAME Atan2
  #define FUNCTION atan2
  #define PRIMAL_FUNCTION atan2
  #include "binaryExpression.tpp"

  /*
   * Implementation for f(a,b) = pow(a,b)
   */
  /**
   * @brief Helper function which checks if the base of the power is negative.
   *
   * Depending of the global option CheckExpressionArguments the function will call
   * CODI_EXCEPTION if a is smaller than zero.
   *
   * @tparam A Type of the first argument.
   * @tparam B Type of the second argument.
   */
  template<typename Real> CODI_INLINE void checkArgumentsPow(const Real& a) {
    if(CheckExpressionArguments) {
      if( TypeTraits<Real>::getBaseValue(a) < 0.0) {
        CODI_EXCEPTION("Negative base for active exponent in pow function. (Value: %0.15e)", TypeTraits<Real>::getBaseValue(a));
      }
    }
  }
  template<typename Real, typename A, typename B> CODI_INLINE Real gradientA_Pow(const A& a, const B& b, const Real& result) {
    CODI_UNUSED(result);
    checkArgumentsPow(a);
    return b * pow(a, b - 1.0);
  }
  template<typename Real, typename A, typename B> CODI_INLINE Real gradientB_Pow(const A& a, const B& b, const Real& result) {
    CODI_UNUSED(b);
    checkArgumentsPow(a);
    if (a > 0.0) {
      return log(a) * result;
    } else {
      return 0.0;
    }
  }
  template<typename Data, typename Real, typename A, typename B> CODI_INLINE void derv11_Pow(Data& data, const A& a, const B& b, const Real& result) {
    checkArgumentsPow(a.getValue());
    a.calcGradient(data, b.getValue() * pow(a.getValue(), b.getValue() - 1.0));
    if (a.getValue() > 0.0) {
      b.calcGradient(data, log(a.getValue()) * result);
    }
  }
  template<typename Data, typename Real, typename A, typename B> CODI_INLINE void derv11M_Pow(Data& data, const A& a, const B& b, const Real& result, const Real& multiplier) {
    checkArgumentsPow(a.getValue());
    a.calcGradient(data, multiplier * b.getValue() * pow(a.getValue(), b.getValue() - 1.0));
    if (a.getValue() > 0.0) {
      b.calcGradient(data, multiplier * log(a.getValue()) * result);
    }
  }
  template<typename Data, typename Real, typename A> CODI_INLINE void derv10_Pow(Data& data, const A& a, const typename TypeTraits<Real>::PassiveReal& b, const Real& result) {
    CODI_UNUSED(result);
    a.calcGradient(data, b * pow(a.getValue(), b - 1.0));
  }
  template<typename Data, typename Real, typename A> CODI_INLINE void derv10M_Pow(Data& data, const A& a, const typename TypeTraits<Real>::PassiveReal& b, const Real& result, const Real& multiplier) {
    CODI_UNUSED(result);
    a.calcGradient(data, multiplier * b * pow(a.getValue(), b - 1.0));
  }
  template<typename Data, typename Real, typename B> CODI_INLINE void derv01_Pow(Data& data, const typename TypeTraits<Real>::PassiveReal& a, const B& b, const Real& result) {
    checkArgumentsPow(a);
    if (a > 0.0) {
      b.calcGradient(data, log(a) * result);
    }
  }
  template<typename Data, typename Real, typename B> CODI_INLINE void derv01M_Pow(Data& data, const typename TypeTraits<Real>::PassiveReal& a, const B& b, const Real& result, const Real& multiplier) {
    checkArgumentsPow(a);
    if (a > 0.0) {
      b.calcGradient(data, multiplier * log(a) * result);
    }
  }
  using std::pow;
  #define NAME Pow
  #define FUNCTION pow
  #define PRIMAL_FUNCTION pow
  #include "binaryExpression.tpp"

  /*
   * Implementation for f(a,b) = Min(a,b)
   */
  template<typename Real, typename A, typename B> CODI_INLINE const typename TypeTraits<Real>::PassiveReal gradientA_Min(const A& a, const B& b, const Real& result) {
    CODI_UNUSED(result);
    if(a < b) {
      return 1.0;
    } else {
      return 0.0;
    }
  }
  template<typename Real, typename A, typename B> CODI_INLINE const typename TypeTraits<Real>::PassiveReal gradientB_Min(const A& a, const B& b, const Real& result) {
    CODI_UNUSED(result);
    if(a < b) {
      return 0.0;
    } else {
      return 1.0;
    }
  }
  template<typename Data, typename Real, typename A, typename B> CODI_INLINE void derv11_Min(Data& data, const A& a, const B& b, const Real& result) {
    CODI_UNUSED(result);
    if(a.getValue() < b.getValue()) {
      a.calcGradient(data);
    } else {
      b.calcGradient(data);
    }
  }
  template<typename Data, typename Real, typename A, typename B> CODI_INLINE void derv11M_Min(Data& data, const A& a, const B& b, const Real& result, const Real& multiplier) {
    CODI_UNUSED(result);
    if(a.getValue() < b.getValue()) {
      a.calcGradient(data, multiplier);
    } else {
      b.calcGradient(data, multiplier);
    }
  }
  template<typename Data, typename Real, typename A> CODI_INLINE void derv10_Min(Data& data, const A& a, const typename TypeTraits<Real>::PassiveReal& b, const Real& result) {
    CODI_UNUSED(result);
    if(a.getValue() < b) {
      a.calcGradient(data);
    }
  }
  template<typename Data, typename Real, typename A> CODI_INLINE void derv10M_Min(Data& data, const A& a, const typename TypeTraits<Real>::PassiveReal& b, const Real& result, const Real& multiplier) {
    CODI_UNUSED(result);
    if(a.getValue() < b) {
      a.calcGradient(data, multiplier);
    }
  }
  template<typename Data, typename Real, typename B> CODI_INLINE void derv01_Min(Data& data, const typename TypeTraits<Real>::PassiveReal& a, const B& b, const Real& result) {
    CODI_UNUSED(result);
    if(a >= b.getValue()) {
      b.calcGradient(data);
    }
  }
  template<typename Data, typename Real, typename B> CODI_INLINE void derv01M_Min(Data& data, const typename TypeTraits<Real>::PassiveReal& a, const B& b, const Real& result, const Real& multiplier) {
    CODI_UNUSED(result);
    if(a >= b.getValue()) {
      b.calcGradient(data, multiplier);
    }
  }
  using std::min;
  #define NAME Min
  #define FUNCTION min
  #define PRIMAL_FUNCTION min
  #include "binaryExpression.tpp"

  /*
   * Forward of fmin to min
   */
  /**
   * @brief Overload for fmin with the CoDiPack expressions.
   *
   * @param[in] a  The first argument of the operation.
   * @param[in] b  The second argument of the operation.
   *
   * @return The implementing expression Min11.
   *
   * @tparam Real  The real type used in the active types.
   * @tparam    A  The expression for the first argument of the function
   * @tparam    B  The expression for the second argument of the function
   */
  template <typename Real, class A, class B>
  CODI_INLINE Min11<Real, A, B> fmin(const Expression<Real, A>& a, const Expression<Real, B>& b) {
    return Min11<Real, A, B>(a.cast(), b.cast());
  }
  /**
   * @brief Overload for fmin with the CoDiPack expressions.
   *
   * @param[in] a  The first argument of the operation.
   * @param[in] b  The second argument of the operation.
   *
   * @return The implementing expression Min10.
   *
   * @tparam Real  The real type used in the active types.
   * @tparam    A  The expression for the first argument of the function
   */
  template <typename Real, class A>
  CODI_INLINE Min10<Real, A> fmin(const Expression<Real, A>& a, const typename TypeTraits<Real>::PassiveReal& b) {
    return Min10<Real, A>(a.cast(), b);
  }
  /**
   * @brief Overload for fmin with the CoDiPack expressions.
   *
   * @param[in] a  The first argument of the operation.
   * @param[in] b  The second argument of the operation.
   *
   * @return The implementing expression Min01.
   *
   * @tparam Real  The real type used in the active types.
   * @tparam    B  The expression for the second argument of the function
   */
  template <typename Real, class B>
  CODI_INLINE Min01<Real, B> fmin(const typename TypeTraits<Real>::PassiveReal& a, const Expression<Real, B>& b) {
    return Min01<Real, B>(a, b.cast());
  }

  /*
   * Implementation for f(a,b) = Max(a,b)
   */
  template<typename Real, typename A, typename B> CODI_INLINE const typename TypeTraits<Real>::PassiveReal gradientA_Max(const A& a, const B& b, const Real& result) {
    CODI_UNUSED(result);
    if(a > b) {
      return 1.0;
    } else {
      return 0.0;
    }
  }
  template<typename Real, typename A, typename B> CODI_INLINE const typename TypeTraits<Real>::PassiveReal gradientB_Max(const A& a, const B& b, const Real& result) {
    CODI_UNUSED(result);
    if(a > b) {
      return 0.0;
    } else {
      return 1.0;
    }
  }
  template<typename Data, typename Real, typename A, typename B> CODI_INLINE void derv11_Max(Data& data, const A& a, const B& b, const Real& result) {
    CODI_UNUSED(result);
    if(a.getValue() > b.getValue()) {
      a.calcGradient(data);
    } else {
      b.calcGradient(data);
    }
  }
  template<typename Data, typename Real, typename A, typename B> CODI_INLINE void derv11M_Max(Data& data, const A& a, const B& b, const Real& result, const Real& multiplier) {
    CODI_UNUSED(result);
    if(a.getValue() > b.getValue()) {
      a.calcGradient(data, multiplier);
    } else {
      b.calcGradient(data, multiplier);
    }
  }
  template<typename Data, typename Real, typename A> CODI_INLINE void derv10_Max(Data& data, const A& a, const typename TypeTraits<Real>::PassiveReal& b, const Real& result) {
    CODI_UNUSED(result);
    if(a.getValue() > b) {
      a.calcGradient(data);
    }
  }
  template<typename Data, typename Real, typename A> CODI_INLINE void derv10M_Max(Data& data, const A& a, const typename TypeTraits<Real>::PassiveReal& b, const Real& result, const Real& multiplier) {
    CODI_UNUSED(result);
    if(a.getValue() > b) {
      a.calcGradient(data, multiplier);
    }
  }
  template<typename Data, typename Real, typename B> CODI_INLINE void derv01_Max(Data& data, const typename TypeTraits<Real>::PassiveReal& a, const B& b, const Real& result) {
    CODI_UNUSED(result);
    if(a <= b.getValue()) {
      b.calcGradient(data);
    }
  }
  template<typename Data, typename Real, typename B> CODI_INLINE void derv01M_Max(Data& data, const typename TypeTraits<Real>::PassiveReal& a, const B& b, const Real& result, const Real& multiplier) {
    CODI_UNUSED(result);
    if(a <= b.getValue()) {
      b.calcGradient(data, multiplier);
    }
  }
  using std::max;
  #define NAME Max
  #define FUNCTION max
  #define PRIMAL_FUNCTION max
  #include "binaryExpression.tpp"

  /*
   * Implementation for f(a,b) = copysign(a,b)
   */
  template<typename Real, typename A, typename B> CODI_INLINE const typename TypeTraits<Real>::PassiveReal gradientA_Copysign(const A& a, const B& b, const Real& result) {
    CODI_UNUSED(result);
    if(a < 0.0) {
      if (b < 0.0) {return (typename TypeTraits<Real>::PassiveReal)1.0;}
      else {return (typename  TypeTraits<Real>::PassiveReal)-1.0;}
    } else if(a > 0.0) {
      if (b < 0.0) {return (typename TypeTraits<Real>::PassiveReal)-1.0;}
      else {return (typename TypeTraits<Real>::PassiveReal)1.0;}
    } else {
      return (typename TypeTraits<Real>::PassiveReal)0.0;
    }
  }

  template<typename Real, typename A, typename B> CODI_INLINE const typename TypeTraits<Real>::PassiveReal gradientB_Copysign(const A& a, const B& b, const Real& result) {
    CODI_UNUSED(a);
    CODI_UNUSED(b);
    CODI_UNUSED(result);
    return 0.0;
  }

  template<typename Data, typename Real, typename A, typename B> CODI_INLINE void derv11_Copysign(Data& data, const A& a, const B& b, const Real& result) {
    CODI_UNUSED(result);
    if(a.getValue() < 0.0) {
      if (b.getValue() < 0.0) {return a.calcGradient(data);}
      else {return a.calcGradient(data, -1.0);}
    } else if(a.getValue() > 0.0) {
      if (b.getValue() < 0.0) {return a.calcGradient(data, -1.0);}
      else {return a.calcGradient(data);}
    } else {
      return a.calcGradient(data, 0.0);
    }
    b.calcGradient(data, 0.0);
  }

  template<typename Data, typename Real, typename A, typename B> CODI_INLINE void derv11M_Copysign(Data& data, const A& a, const B& b, const Real& result, const Real& multiplier) {
    CODI_UNUSED(result);
    if(a.getValue() < 0.0) {
      if (b.getValue() < 0.0) {return a.calcGradient(data, multiplier);}
      else {return a.calcGradient(data, -1.0*multiplier);}
    } else if(a.getValue() > 0.0) {
      if (b.getValue() < 0.0) {return a.calcGradient(data, -1.0*multiplier);}
      else {return a.calcGradient(data, multiplier);}
    } else {
      return a.calcGradient(data, 0.0);
    }
    b.calcGradient(data, 0.0);
  }

  template<typename Data, typename Real, typename A> CODI_INLINE void derv10_Copysign(Data& data, const A& a, const typename TypeTraits<Real>::PassiveReal& b, const Real& result) {
    CODI_UNUSED(result);
    if(a.getValue() < 0.0) {
      if (b < 0.0) {return a.calcGradient(data);}
      else {return a.calcGradient(data, -1.0);}
    } else if(a.getValue() > 0.0) {
      if (b < 0.0) {return a.calcGradient(data, -1.0);}
      else {return a.calcGradient(data);}
    } else {
      return a.calcGradient(data, 0.0);
    }
  }

  template<typename Data, typename Real, typename A> CODI_INLINE void derv10M_Copysign(Data& data, const A& a, const typename TypeTraits<Real>::PassiveReal& b, const Real& result, const Real& multiplier) {
    CODI_UNUSED(result);
    if(a.getValue() < 0.0) {
      if (b < 0.0) {return a.calcGradient(data, multiplier);}
      else {return a.calcGradient(data, -1.0*multiplier);}
    } else if(a.getValue() > 0.0) {
      if (b < 0.0) {return a.calcGradient(data, -1.0*multiplier);}
      else {return a.calcGradient(data, multiplier);}
    } else {
      return a.calcGradient(data, 0.0);
    }
  }

  template<typename Data, typename Real, typename B> CODI_INLINE void derv01_Copysign(Data& data, const typename TypeTraits<Real>::PassiveReal& a, const B& b, const Real& result) {
    CODI_UNUSED(a);
    CODI_UNUSED(b);
    CODI_UNUSED(result);
    b.calcGradient(data, 0.0);
  }

  template<typename Data, typename Real, typename B> CODI_INLINE void derv01M_Copysign(Data& data, const typename TypeTraits<Real>::PassiveReal& a, const B& b, const Real& result, const Real& multiplier) {
    CODI_UNUSED(a);
    CODI_UNUSED(b);
    CODI_UNUSED(result);
    b.calcGradient(data, 0.0);
  }

  using std::copysign;
  #define NAME Copysign
  #define FUNCTION copysign
  #define PRIMAL_FUNCTION copysign
  #include "binaryExpression.tpp"

  /**
   * @brief Overload for copysignf with the CoDiPack expressions.
   *
   * @param[in] a  The first argument of the operation.
   * @param[in] b  The second argument of the operation.
   *
   * @return The implementing expression Copysign11.
   *
   * @tparam Real  The real type used in the active types.
   * @tparam    A  The expression for the first argument of the function
   * @tparam    B  The expression for the second argument of the function
   */
  template <typename Real, class A, class B>
  CODI_INLINE Copysign11<Real, A, B> copysignf(const Expression<Real, A>& a, const Expression<Real, B>& b) {
    return Copysign11<Real, A, B>(a.cast(), b.cast());
  }

  /**
   * @brief Overload for copysignf with the CoDiPack expressions.
   *
   * @param[in] a  The first argument of the operation.
   * @param[in] b  The second argument of the operation.
   *
   * @return The implementing expression Copysign10.
   *
   * @tparam Real  The real type used in the active types.
   * @tparam    A  The expression for the first argument of the function
   * @tparam    B  The expression for the second argument of the function
   */
  template <typename Real, class A>
  CODI_INLINE Copysign10<Real, A> copysignf(const Expression<Real, A>& a, const typename TypeTraits<Real>::PassiveReal& b) {
    return Copysign10<Real, A>(a.cast(), b);
  }

  /**
   * @brief Overload for copysignf with the CoDiPack expressions.
   *
   * @param[in] a  The first argument of the operation.
   * @param[in] b  The second argument of the operation.
   *
   * @return The implementing expression Copysign01.
   *
   * @tparam Real  The real type used in the active types.
   * @tparam    A  The expression for the first argument of the function
   * @tparam    B  The expression for the second argument of the function
   */
  template <typename Real, class B>
  CODI_INLINE Copysign01<Real, B> copysignf(const typename TypeTraits<Real>::PassiveReal& a, const Expression<Real, B>& b) {
    return Copysign01<Real, B>(a, b.cast());
  }

  /*
   * Forward of fmax to max
   */
  /**
   * @brief Overload for fmax with the CoDiPack expressions.
   *
   * @param[in] a  The first argument of the operation.
   * @param[in] b  The second argument of the operation.
   *
   * @return The implementing expression Max11.
   *
   * @tparam Real  The real type used in the active types.
   * @tparam    A  The expression for the first argument of the function
   * @tparam    B  The expression for the second argument of the function
   */
  template <typename Real, class A, class B>
  CODI_INLINE Max11<Real, A, B> fmax(const Expression<Real, A>& a, const Expression<Real, B>& b) {
    return Max11<Real, A, B>(a.cast(), b.cast());
  }

  /**
   * @brief Overload for fmax with the CoDiPack expressions.
   *
   * @param[in] a  The first argument of the operation.
   * @param[in] b  The second argument of the operation.
   *
   * @return The implementing expression Max10.
   *
   * @tparam Real  The real type used in the active types.
   * @tparam    A  The expression for the first argument of the function
   */
  template <typename Real, class A>
  CODI_INLINE Max10<Real, A> fmax(const Expression<Real, A>& a, const typename TypeTraits<Real>::PassiveReal& b) {
    return Max10<Real, A>(a.cast(), b);
  }
  /**
   * @brief Overload for fmax with the CoDiPack expressions.
   *
   * @param[in] a  The first argument of the operation.
   * @param[in] b  The second argument of the operation.
   *
   * @return The implementing expression Max01.
   *
   * @tparam Real  The real type used in the active types.
   * @tparam    B  The expression for the second argument of the function
   */
  template <typename Real, class B>
  CODI_INLINE Max01<Real, B> fmax(const typename TypeTraits<Real>::PassiveReal& a, const Expression<Real, B>& b) {
    return Max01<Real, B>(a, b.cast());
  }


  #undef CODI_OPERATOR_HELPER

  /*
   * Conditional operators should behave exactly the same as with
   * non-active arguments so in each of the cases below the getValue()
   * function is called to extract the value of the expression
   */
  #define CODI_DEFINE_CONDITIONAL(OPERATOR, OP) \
    /** @brief Overload for OP with the CoDiPack expressions. @param[in] a The first argument of the operation. @param[in] b The second argument of the operation. @return The operation returns the same value the same version with double arguments. @tparam Real The real type used in the active types. @tparam A The expression for the first argument of the function @tparam B The expression for the second argument of the function*/ \
    template<typename Real, class A, class B> \
    CODI_INLINE bool OPERATOR(const Expression<Real, A>& a, const Expression<Real, B>& b) { \
      return a.getValue() OP b.getValue(); \
    } \
    \
    /** @brief Overload for OP with the CoDiPack expressions. @param[in] a The first argument of the operation. @param[in] b The second argument of the operation. @return The operation returns the same value the same version with double arguments. @tparam Real The real type used in the active types. @tparam A The expression for the first argument of the function */ \
    template<typename Real, class A> \
    CODI_INLINE bool OPERATOR(const Expression<Real, A>& a, const typename TypeTraits<Real>::PassiveReal& b) { \
      return a.getValue() OP b; \
    } \
    \
    /** @brief Overload for OP with the CoDiPack expressions. @param[in] a The first argument of the operation. @param[in] b The second argument of the operation. @return The operation returns the same value the same version with double arguments. @tparam Real The real type used in the active types. @tparam B The expression for the second argument of the function*/ \
    template<typename Real, class B> \
    CODI_INLINE bool OPERATOR(const typename TypeTraits<Real>::PassiveReal& a, const Expression<Real, B>& b) { \
      return a OP b.getValue(); \
    } \
    /** @brief Overload for OP with the CoDiPack expressions. @param[in] a The first argument of the operation. @param[in] b The second argument of the operation. @return The operation returns the same value the same version with double arguments. @tparam Real The real type used in the active types. @tparam A The expression for the first argument of the function */ \
    template<typename Real, class A> \
    CODI_INLINE bool OPERATOR(const Expression<Real, A>& a, const int& b) { \
      return a.getValue() OP b; \
    } \
    \
    /** @brief Overload for OP with the CoDiPack expressions. @param[in] a The first argument of the operation. @param[in] b The second argument of the operation. @return The operation returns the same value the same version with double arguments. @tparam Real The real type used in the active types. @tparam B The expression for the second argument of the function*/ \
    template<typename Real, class B>            \
    CODI_INLINE bool OPERATOR(const int& a, const Expression<Real, B>& b) { \
      return a OP b.getValue(); \
    } \
    \
    /** @brief Overload for OP with the CoDiPack expressions. @param[in] a The first argument of the operation. @param[in] b The second argument of the operation. @return The operation returns the same value the same version with double arguments. @tparam Real The real type used in the active types. @tparam A The expression for the first argument of the function */ \
    template<typename Real, class A> \
    CODI_INLINE bool OPERATOR(const Expression<Real, A>& a, const long& b) { \
      return a.getValue() OP b; \
    } \
    \
    /** @brief Overload for OP with the CoDiPack expressions. @param[in] a The first argument of the operation. @param[in] b The second argument of the operation. @return The operation returns the same value the same version with double arguments. @tparam Real The real type used in the active types. @tparam B The expression for the second argument of the function*/ \
    template<typename Real, class B>            \
    CODI_INLINE bool OPERATOR(const long& a, const Expression<Real, B>& b) { \
      return a OP b.getValue(); \
    }

  CODI_DEFINE_CONDITIONAL(operator==, ==)
  CODI_DEFINE_CONDITIONAL(operator!=, !=)
  CODI_DEFINE_CONDITIONAL(operator>, >)
  CODI_DEFINE_CONDITIONAL(operator<, <)
  CODI_DEFINE_CONDITIONAL(operator>=, >=)
  CODI_DEFINE_CONDITIONAL(operator<=, <=)
  CODI_DEFINE_CONDITIONAL(operator&&, &&)
  CODI_DEFINE_CONDITIONAL(operator||, ||)

  #undef CODI_DEFINE_CONDITIONAL

  #define CODI_DEFINE_UNARY_CONDITIONAL(OPERATOR, OP) \
    /** @brief Overload for OP with the CoDiPack expressions. @param[in] a The argument of the operation. @return The operation returns the same value the same version with double arguments. @tparam Real The real type used in the active types. @tparam A The expression for the first argument of the function */ \
    template<typename Real, class A> \
    CODI_INLINE bool OPERATOR(const Expression<Real, A>& a) { \
      return OP a.getValue(); \
    }

  CODI_DEFINE_UNARY_CONDITIONAL(operator!, !)

  #undef CODI_DEFINE_UNARY_CONDITIONAL

  #define CODI_OPERATOR_HELPER(NAME, OP) \
    /** @brief Helper function to call operators as a function. @param[in] a The argument of the operation. @return The value of OP b @tparam A The expression for the argument of the function*/ \
    template<typename A> \
    CODI_INLINE auto primal_ ## NAME(const A& a) -> decltype(OP a) { \
      return OP a; \
    }

  /*
   * Now all unary functions are implemented.
   *
   * We use the naming scheme gradName for the gradient computation of the function.
   *
   */

  template<typename Real> CODI_INLINE typename TypeTraits<Real>::PassiveReal gradient_UnaryMinus(const Real& a, const Real& result) {
    CODI_UNUSED(a);
    CODI_UNUSED(result);
    return -1.0;
  }
  CODI_OPERATOR_HELPER(UnaryMinus, -)
  #define NAME UnaryMinus
  #define FUNCTION operator -
  #define PRIMAL_FUNCTION primal_UnaryMinus
  #include "unaryExpression.tpp"

  template<typename Real> CODI_INLINE Real gradient_Sqrt(const Real& a, const Real& result) {
    if(CheckExpressionArguments) {
      if(0.0 > TypeTraits<Real>::getBaseValue(a)) {
        CODI_EXCEPTION("Sqrt of negative value or zero.(Value: %0.15e)", TypeTraits<Real>::getBaseValue(a));
      }
    }
    if(result != 0.0) {
      return 0.5 / result;
    } else {
      return (Real)0.0;
    }
  }
  using std::sqrt;
  #define NAME Sqrt
  #define FUNCTION sqrt
  #define PRIMAL_FUNCTION sqrt
  #include "unaryExpression.tpp"

  template<typename Real> CODI_INLINE Real gradient_Cbrt(const Real& a, const Real& result) {
    if(CheckExpressionArguments) {
      if(0.0 == TypeTraits<Real>::getBaseValue(a)) {
        CODI_EXCEPTION("Cbrt of zero value.(Value: %0.15e)", TypeTraits<Real>::getBaseValue(a));
      }
    }
    if(result != 0.0) {
      return 1.0 / (3.0 * result * result);
    } else {
      return (Real)0.0;
    }
  }
  using std::cbrt;
  #define NAME Cbrt
  #define FUNCTION cbrt
  #define PRIMAL_FUNCTION cbrt
  #include "unaryExpression.tpp"

  template<typename Real> CODI_INLINE Real gradient_Tanh(const Real& a, const Real& result) {
    CODI_UNUSED(a);
    return 1 - result * result;
  }
  using std::tanh;
  #define NAME Tanh
  #define FUNCTION tanh
  #define PRIMAL_FUNCTION tanh
  #include "unaryExpression.tpp"

  template<typename Real> CODI_INLINE Real gradient_Log(const Real& a, const Real& result) {
    CODI_UNUSED(result);
    if(CheckExpressionArguments) {
      if(0.0 > TypeTraits<Real>::getBaseValue(a)) {
        CODI_EXCEPTION("Logarithm of negative value or zero.(Value: %0.15e)", TypeTraits<Real>::getBaseValue(a));
      }
    }
    return 1.0 / a;
  }
  using std::log;
  #define NAME Log
  #define FUNCTION log
  #define PRIMAL_FUNCTION log
  #include "unaryExpression.tpp"

  template<typename Real> CODI_INLINE Real gradient_Log10(const Real& a, const Real& result) {
    CODI_UNUSED(result);
    if(CheckExpressionArguments) {
      if(0.0 > TypeTraits<Real>::getBaseValue(a)) {
        CODI_EXCEPTION("Logarithm of negative value or zero.(Value: %0.15e)", TypeTraits<Real>::getBaseValue(a));
      }
    }
    return 0.434294481903252 / a;
  }
  using std::log10;
  #define NAME Log10
  #define FUNCTION log10
  #define PRIMAL_FUNCTION log10
  #include "unaryExpression.tpp"

  template<typename Real> CODI_INLINE Real gradient_Sin(const Real& a, const Real& result) {
    CODI_UNUSED(result);
    return cos(a);
  }
  using std::sin;
  #define NAME Sin
  #define FUNCTION sin
  #define PRIMAL_FUNCTION sin
  #include "unaryExpression.tpp"

  template<typename Real> CODI_INLINE Real gradient_Cos(const Real& a, const Real& result) {
    CODI_UNUSED(result);
    return -sin(a);
  }
  using std::cos;
  #define NAME Cos
  #define FUNCTION cos
  #define PRIMAL_FUNCTION cos
  #include "unaryExpression.tpp"

  using std::sqrt;
  template<typename Real> CODI_INLINE Real gradient_Asin(const Real& a, const Real& result) {
    CODI_UNUSED(result);
    if(CheckExpressionArguments) {
      if(TypeTraits<Real>::getBaseValue(a) <= -1.0 || 1.0 <= TypeTraits<Real>::getBaseValue(a)) {
        CODI_EXCEPTION("asin outside of (-1, 1).(Value: %0.15e)", TypeTraits<Real>::getBaseValue(a));
      }
    }
    return 1.0 / sqrt(1.0 - a * a);
  }
  using std::asin;
  #define NAME Asin
  #define FUNCTION asin
  #define PRIMAL_FUNCTION asin
  #include "unaryExpression.tpp"

  template<typename Real> CODI_INLINE Real gradient_Acos(const Real& a, const Real& result) {
    CODI_UNUSED(result);
    if(CheckExpressionArguments) {
      if(TypeTraits<Real>::getBaseValue(a) <= -1.0 || 1.0 <= TypeTraits<Real>::getBaseValue(a)) {
        CODI_EXCEPTION("acos outside of (-1, 1).(Value: %0.15e)", TypeTraits<Real>::getBaseValue(a));
      }
    }
    return -1.0 / sqrt(1.0 - a * a);
  }
  using std::acos;
  #define NAME Acos
  #define FUNCTION acos
  #define PRIMAL_FUNCTION acos
  #include "unaryExpression.tpp"

  template<typename Real> CODI_INLINE Real gradient_Atan(const Real& a, const Real& result) {
    CODI_UNUSED(result);
    return 1.0 / (1 + a * a);
  }
  using std::atan;
  #define NAME Atan
  #define FUNCTION atan
  #define PRIMAL_FUNCTION atan
  #include "unaryExpression.tpp"

  template<typename Real> CODI_INLINE Real gradient_Sinh(const Real& a, const Real& result) {
    CODI_UNUSED(result);
    return cosh(a);
  }
  using std::sinh;
  #define NAME Sinh
  #define FUNCTION sinh
  #define PRIMAL_FUNCTION sinh
  #include "unaryExpression.tpp"

  template<typename Real> CODI_INLINE Real gradient_Cosh(const Real& a, const Real& result) {
    CODI_UNUSED(result);
    return sinh(a);
  }
  using std::cosh;
  #define NAME Cosh
  #define FUNCTION cosh
  #define PRIMAL_FUNCTION cosh
  #include "unaryExpression.tpp"

  template<typename Real> CODI_INLINE Real gradient_Exp(const Real& a, const Real& result) {
    CODI_UNUSED(a);
    return result;
  }
  using std::exp;
  #define NAME Exp
  #define FUNCTION exp
  #define PRIMAL_FUNCTION exp
  #include "unaryExpression.tpp"

  template<typename Real> CODI_INLINE Real gradient_Atanh(const Real& a, const Real& result) {
    CODI_UNUSED(result);
    if(CheckExpressionArguments) {
      if(TypeTraits<Real>::getBaseValue(a) <= -1.0 || 1.0 <= TypeTraits<Real>::getBaseValue(a)) {
        CODI_EXCEPTION("atanh outside of (-1, 1).(Value: %0.15e)", TypeTraits<Real>::getBaseValue(a));
      }
    }
    return 1.0 / (1 - a * a);
  }
  using std::atanh;
  #define NAME Atanh
  #define FUNCTION atanh
  #define PRIMAL_FUNCTION atanh
  #include "unaryExpression.tpp"

  template<typename Real> CODI_INLINE Real gradient_Abs(const Real& a, const Real& result) {
    CODI_UNUSED(result);
    if(a < 0.0) {
      return (Real)-1.0;
    } else if(a > 0.0) {
      return (Real)1.0;
    } else {
      return (Real)0.0;
    }
  }
  using std::abs;
  #define NAME Abs
  #define FUNCTION abs
  #define PRIMAL_FUNCTION abs
  #include "unaryExpression.tpp"

  template<typename Real> CODI_INLINE Real gradient_Tan(const Real& a, const Real& result) {
    CODI_UNUSED(result);
    if(CheckExpressionArguments) {
      if(0.0 == cos(TypeTraits<Real>::getBaseValue(a))) {
        CODI_EXCEPTION("Tan evaluated at (0.5  + i) * PI.(Value: %0.15e)", TypeTraits<Real>::getBaseValue(a));
      }
    }
    Real tmp = 1.0 / cos(a);
    return tmp * tmp;
  }
  using std::tan ;
  #define NAME Tan
  #define FUNCTION tan
  #define PRIMAL_FUNCTION tan
  #include "unaryExpression.tpp"

  template<typename Real> CODI_INLINE Real gradient_Erf(const Real& a, const Real& result) {
    CODI_UNUSED(result);
    return 1.128379167095513 * exp( -(a * a) ); // erf'(a) = 2.0 / sqrt(pi) * exp(-a^2)
  }
  using std::erf;
  #define NAME Erf
  #define FUNCTION erf
  #define PRIMAL_FUNCTION erf
  #include "unaryExpression.tpp"

  template<typename Real> CODI_INLINE Real gradient_Erfc(const Real& a, const Real& result) {
    CODI_UNUSED(result);
    return -1.128379167095513 * exp( -(a * a) ); // erfc'(a) = - 2.0 / sqrt(pi) * exp(-a^2)
  }
  using std::erfc;
  #define NAME Erfc
  #define FUNCTION erfc
  #define PRIMAL_FUNCTION erfc
  #include "unaryExpression.tpp"

  template<typename Real> CODI_INLINE Real gradient_Tgamma(const Real& a, const Real& result) {
    if(a <= 0.0) {
      std::cout << "Derivative for gamma function only for positive arguments at the moment" << std::endl;
      std::exit(1);
    }

    // Implementation of the digamma function is taken from John Burkardt,
    // http://people.sc.fsu.edu/~jburkardt/cpp_src/asa103/asa103.cpp
    //
    // Definition of Gamma(a): https://en.wikipedia.org/wiki/Gamma_function
    // Definition of DiGamma(a): https://en.wikipedia.org/wiki/Digamma_function
    // Differentation is Gamma'(a) = Gamma(a) * DiGamma(a)

    Real diGamma = 0.0;
    if(a <= 0.000001) { // special case for small numbers
      const Real eulerMascheroni = 0.57721566490153286060;
      diGamma = -eulerMascheroni - 1.0/a + 1.6449340668482264365*a;
    } else {
      // shift DiGamma(a) = DiGamma(a + 1) - 1/a
      // we require a large such that the approximation below is more accurate
      Real shiftBound = 8.5;

      Real shiftedValue = a;
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
  using std::tgamma;
  #define NAME Tgamma
  #define FUNCTION tgamma
  #define PRIMAL_FUNCTION tgamma
  #include "unaryExpression.tpp"

  #undef CODI_OPERATOR_HELPER

  /**
   * @brief The fabs function is redirected to abs.
   *
   * @param[in] a The argument of the operation.
   *
   * @return The implementation of Abs.
   *
   * @tparam Real The real type used in the active types.
   * @tparam A The expression for the argument of the function
   */
  template<typename Real, class A>
  CODI_INLINE Abs<Real, A> fabs(const codi::Expression<Real, A>& a) {
    return Abs<Real, A>(a.cast());
  }

  /**
   * @brief Overload for operator+ with the CoDiPack expressions.
   *
   * @param[in] a The argument of the expression.
   * @return The argument a.
   *
   * @tparam Real The real type used in the active types.
   * @tparam A The expression for the argument of the function.
   */
  template<typename Real, class A>
  CODI_INLINE const Expression<Real, A>& operator+(const Expression<Real, A>& a) {
    return a;
  }

  /***************************************************************************************
   * Functions that do not need derivatives.
   ****************************************************************************************/
  using std::isinf;
  /**
   * @brief Overload for the isinf function with expressions.
   *
   * @param[in] a The argument of the function.
   *
   * @return The result of isinf on the primal value.
   *
   * @tparam Real The real type used in the active types.
   * @tparam A The expression for the argument of the function
   */
  template<typename Real, class A>
  CODI_INLINE bool isinf(const codi::Expression<Real, A>& a) {
    return isinf(a.getValue());
  }

  using std::isnan;
  /**
   * @brief Overload for the isnan function with expressions.
   *
   * @param[in] a The argument of the function.
   *
   * @return The result of isnan on the primal value.
   *
   * @tparam Real The real type used in the active types.
   * @tparam A The expression for the argument of the function
   */
  template<typename Real, class A>
  CODI_INLINE bool isnan(const codi::Expression<Real, A>& a) {
    return isnan(a.getValue());
  }

  using std::floor;
  /**
   * @brief Overload for the floor function with expressions.
   *
   * @param[in] a The argument of the function.
   *
   * @return The result of floor on the primal value.
   *
   * @tparam Real The real type used in the active types.
   * @tparam A The expression for the argument of the function
   */
  template<typename Real, class A>
  CODI_INLINE typename codi::TypeTraits<Real>::PassiveReal floor(const codi::Expression<Real, A>& a) {
    return floor(a.getValue());
  }

  using std::ceil;
  /**
   * @brief Overload for the floor function with expressions.
   *
   * @param[in] a The argument of the function.
   *
   * @return The result of floor on the primal value.
   *
   * @tparam Real The real type used in the active types.
   * @tparam A The expression for the argument of the function
   */
  template<typename Real, class A>
  CODI_INLINE typename codi::TypeTraits<Real>::PassiveReal ceil(const codi::Expression<Real, A>& a) {
    return ceil(a.getValue());
  }

}
