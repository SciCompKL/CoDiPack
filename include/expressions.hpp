/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015 Chair for Scientific Computing (SciComp), TU Kaiserslautern
 * Homepage: http://www.scicomp.uni-kl.de
 * Contact:  Prof. Nicolas R. Gauger (codi@scicomp.uni-kl.de)
 *
 * Lead developers: Max Sagebaum, Tim Albring (SciComp, TU Kaiserslautern)
 *
 * This file is part of CoDiPack (http://www.scicomp.uni-kl.de/software/codi).
 *
 * CoDiPack is free software: you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation, either version 2 of the
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
    inline const A& cast() const {
      return static_cast<const A&>(*this);
    }

    /**
     * @brief Calculate the gradient of the expression.
     *
     * Calculate the gradient of the mathematical operation that this
     * expression represents and pass the result to its argument.
     * For functions f(a), pass df/da to the argument.
     *
     * @param[inout] data A helper value which the tape can define and use for the evaluation.
     *
     * @tparam Data The type for the tape data.
     */
    template<typename Data>
    inline void calcGradient(Data& data) const {
      cast().calcGradient(data);
    }

    /**
     * @brief Calculate the gradient of the expression.
     *
     * Calculate the gradient of the mathematical operation that this
     * expression represents and pass the result to its argument.
     * For functions f(a), pass multiplier * df/da to the argument.
     *
     * @param[inout]     data A helper value which the tape can define and use for the evaluation.
     * @param[in]  multiplier The Jacobi from the expression where this expression was used as an argument.
     *
     * @tparam Data The type for the tape data.
     */
    template<typename Data>
    inline void calcGradient(Data& data, const Real& multiplier) const {
      cast().calcGradient(data, multiplier);
    }

    /**
     * @brief Return the numerical value of the expression.
     *
     * @return The value of the expression.
     */
    inline const Real getValue() const {
      return cast().getValue();
    }

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
   * inline auto primal_NAME(const A& a, const B& b);
   *
   * @param   NAME  The name for the generated function.
   * @param   OP    The sign of the operator that the function calls.
   */
  #define CODI_OPERATOR_HELPER(NAME, OP) \
    /** @brief Helper function to call operators as a function. @param[in] a The first argument of the operation. @param[in] b The second argument of the operation. @return The value of a OP b @tparam A The expression for the first argument of the function @tparam B The expression for the second argument of the function*/ \
    template<typename A, typename B> \
    inline auto primal_ ## NAME(const A& a, const B& b) -> decltype(a OP b) { \
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
   * There is no implementation for 00 because no variable is active and thus the derivative would be zero.
   *
   * If the M is present the method is implemented with the multiplier as an argument.
   */

  /*
   * Implementation for f(a,b) = a + b
   * df/da = 1 and likewise for df/db so simply
   * call a and b's versions of calcGradient
   */
  template<typename Data, typename Real, typename A, typename B> inline void derv11_Add(Data& data, const A& a, const B& b, const Real& result) {
    CODI_UNUSED(result);
    a.calcGradient(data);
    b.calcGradient(data);
  }
  template<typename Data, typename Real, typename A, typename B> inline void derv11M_Add(Data& data, const A& a, const B& b, const Real& result, const Real& multiplier) {
    a.calcGradient(data, multiplier);
    b.calcGradient(data, multiplier);
  }
  template<typename Data, typename Real, typename A> inline void derv10_Add(Data& data, const A& a, const typename TypeTraits<Real>::PassiveReal& b, const Real& result) {
    CODI_UNUSED(result);
    a.calcGradient(data);
  }
  template<typename Data, typename Real, typename A> inline void derv10M_Add(Data& data, const A& a, const typename TypeTraits<Real>::PassiveReal& b, const Real& result, const Real& multiplier) {
    CODI_UNUSED(result);
    a.calcGradient(data, multiplier);
  }
  template<typename Data, typename Real, typename B> inline void derv01_Add(Data& data, const typename TypeTraits<Real>::PassiveReal& a, const B& b, const Real& result) {
    CODI_UNUSED(result);
    b.calcGradient(data);
  }
  template<typename Data, typename Real, typename B> inline void derv01M_Add(Data& data, const typename TypeTraits<Real>::PassiveReal& a, const B& b, const Real& result, const Real& multiplier) {
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
  template<typename Data, typename Real, typename A, typename B> inline void derv11_Subtract(Data& data, const A& a, const B& b, const Real& result) {
    CODI_UNUSED(result);
    a.calcGradient(data);
    b.calcGradient(data, -1.0);
  }
  template<typename Data, typename Real, typename A, typename B> inline void derv11M_Subtract(Data& data, const A& a, const B& b, const Real& result, const Real& multiplier) {
    CODI_UNUSED(result);
    a.calcGradient(data, multiplier);
    b.calcGradient(data, -multiplier);
  }
  template<typename Data, typename Real, typename A> inline void derv10_Subtract(Data& data, const A& a, const typename TypeTraits<Real>::PassiveReal& b, const Real& result) {
    CODI_UNUSED(result);
    a.calcGradient(data);
  }
  template<typename Data, typename Real, typename A> inline void derv10M_Subtract(Data& data, const A& a, const typename TypeTraits<Real>::PassiveReal& b, const Real& result, const Real& multiplier) {
    CODI_UNUSED(result);
    a.calcGradient(data, multiplier);
  }
  template<typename Data, typename Real, typename B> inline void derv01_Subtract(Data& data, const typename TypeTraits<Real>::PassiveReal& a, const B& b, const Real& result) {
    CODI_UNUSED(result);
    b.calcGradient(data, -1.0);
  }
  template<typename Data, typename Real, typename B> inline void derv01M_Subtract(Data& data, const typename TypeTraits<Real>::PassiveReal& a, const B& b, const Real& result, const Real& multiplier) {
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
  template<typename Data, typename Real, typename A, typename B> inline void derv11_Multiply(Data& data, const A& a, const B& b, const Real& result) {
    CODI_UNUSED(result);
    a.calcGradient(data, b.getValue());
    b.calcGradient(data, a.getValue());
  }
  template<typename Data, typename Real, typename A, typename B> inline void derv11M_Multiply(Data& data, const A& a, const B& b, const Real& result, const Real& multiplier) {
    CODI_UNUSED(result);
    a.calcGradient(data, b.getValue() * multiplier);
    b.calcGradient(data, a.getValue() * multiplier);
  }
  template<typename Data, typename Real, typename A> inline void derv10_Multiply(Data& data, const A& a, const typename TypeTraits<Real>::PassiveReal& b, const Real& result) {
    CODI_UNUSED(result);
    a.calcGradient(data, b);
  }
  template<typename Data, typename Real, typename A> inline void derv10M_Multiply(Data& data, const A& a, const typename TypeTraits<Real>::PassiveReal& b, const Real& result, const Real& multiplier) {
    CODI_UNUSED(result);
    a.calcGradient(data, b * multiplier);
  }
  template<typename Data, typename Real, typename B> inline void derv01_Multiply(Data& data, const typename TypeTraits<Real>::PassiveReal& a, const B& b, const Real& result) {
    CODI_UNUSED(result);
    b.calcGradient(data, a);
  }
  template<typename Data, typename Real, typename B> inline void derv01M_Multiply(Data& data, const typename TypeTraits<Real>::PassiveReal& a, const B& b, const Real& result, const Real& multiplier) {
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
  template<typename Real> inline void checkArgumentsDivide(const Real& b) {
    if(CheckExpressionArguments) {
      if( 0.0 == TypeTraits<Real>::getBaseValue(b)) {
        CODI_EXCEPTION("Division called with divisor of zero.");
      }
    }
  }
  template<typename Data, typename Real, typename A, typename B> inline void derv11_Divide(Data& data, const A& a, const B& b, const Real& result) {
    checkArgumentsDivide(b.getValue());
    Real one_over_b = 1.0 / b.getValue();
    a.calcGradient(data, one_over_b);
    b.calcGradient(data, -result * one_over_b);
  }
  template<typename Data, typename Real, typename A, typename B> inline void derv11M_Divide(Data& data, const A& a, const B& b, const Real& result, const Real& multiplier) {
    checkArgumentsDivide(b.getValue());
    Real one_over_b = multiplier / b.getValue();
    a.calcGradient(data, one_over_b);
    b.calcGradient(data, -result * one_over_b);
  }
  template<typename Data, typename Real, typename A> inline void derv10_Divide(Data& data, const A& a, const typename TypeTraits<Real>::PassiveReal& b, const Real& result) {
    checkArgumentsDivide(b);
    Real one_over_b = 1.0 / b;
    a.calcGradient(data, one_over_b);
  }
  template<typename Data, typename Real, typename A> inline void derv10M_Divide(Data& data, const A& a, const typename TypeTraits<Real>::PassiveReal& b, const Real& result, const Real& multiplier) {
    checkArgumentsDivide(b);
    Real one_over_b = multiplier / b;
    a.calcGradient(data, one_over_b);
  }
  template<typename Data, typename Real, typename B> inline void derv01_Divide(Data& data, const typename TypeTraits<Real>::PassiveReal& a, const B& b, const Real& result) {
    checkArgumentsDivide(b.getValue());
    Real one_over_b = 1.0 / b.getValue();
    b.calcGradient(data, -result * one_over_b);
  }
  template<typename Data, typename Real, typename B> inline void derv01M_Divide(Data& data, const typename TypeTraits<Real>::PassiveReal& a, const B& b, const Real& result, const Real& multiplier) {
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
  template<typename A, typename B> inline void checkArgumentsAtan2(const A& a, const B& b) {
    if(CheckExpressionArguments) {
      if( 0.0 == TypeTraits<A>::getBaseValue(a) &&
          0.0 == TypeTraits<B>::getBaseValue(b)) {
        CODI_EXCEPTION("atan2 called at point (0,0).");
      }
    }
  }
  template<typename Data, typename Real, typename A, typename B> inline void derv11_Atan2(Data& data, const A& a, const B& b, const Real& result) {
    CODI_UNUSED(result);
    checkArgumentsAtan2(a.getValue(), b.getValue());
    Real divisor = a.getValue() * a.getValue() + b.getValue() * b.getValue();
    divisor = 1.0 / divisor;
    a.calcGradient(data, b.getValue() * divisor);
    b.calcGradient(data, -a.getValue() * divisor);
  }
  template<typename Data, typename Real, typename A, typename B> inline void derv11M_Atan2(Data& data, const A& a, const B& b, const Real& result, const Real& multiplier) {
    CODI_UNUSED(result);
    checkArgumentsAtan2(a.getValue(), b.getValue());
    Real divisor = a.getValue() * a.getValue() + b.getValue() * b.getValue();
    divisor = 1.0 / divisor;
    a.calcGradient(data, multiplier * b.getValue() * divisor);
    b.calcGradient(data, multiplier * -a.getValue() * divisor);
  }
  template<typename Data, typename Real, typename A> inline void derv10_Atan2(Data& data, const A& a, const typename TypeTraits<Real>::PassiveReal& b, const Real& result) {
    CODI_UNUSED(result);
    checkArgumentsAtan2(a.getValue(), b);
    Real divisor = a.getValue() * a.getValue() + b * b;
    divisor = 1.0 / divisor;
    a.calcGradient(data, b * divisor);
  }
  template<typename Data, typename Real, typename A> inline void derv10M_Atan2(Data& data, const A& a, const typename TypeTraits<Real>::PassiveReal& b, const Real& result, const Real& multiplier) {
    CODI_UNUSED(result);
    checkArgumentsAtan2(a.getValue(), b);
    Real divisor = a.getValue() * a.getValue() + b * b;
    divisor = 1.0 / divisor;
    a.calcGradient(data, multiplier * b * divisor);
  }
  template<typename Data, typename Real, typename B> inline void derv01_Atan2(Data& data, const typename TypeTraits<Real>::PassiveReal& a, const B& b, const Real& result) {
    CODI_UNUSED(result);
    checkArgumentsAtan2(a, b.getValue());
    Real divisor = a * a + b.getValue() * b.getValue();
    divisor = 1.0 / divisor;
    b.calcGradient(data, -a * divisor);
  }
  template<typename Data, typename Real, typename B> inline void derv01M_Atan2(Data& data, const typename TypeTraits<Real>::PassiveReal& a, const B& b, const Real& result, const Real& multiplier) {
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
  template<typename Real> inline void checkArgumentsPow(const Real& a) {
    if(CheckExpressionArguments) {
      if( TypeTraits<Real>::getBaseValue(a) < 0.0) {
        CODI_EXCEPTION("Negative base for active exponent in pow function. (Value: %0.15e)", TypeTraits<Real>::getBaseValue(a));
      }
    }
  }
  template<typename Data, typename Real, typename A, typename B> inline void derv11_Pow(Data& data, const A& a, const B& b, const Real& result) {
    checkArgumentsPow(a.getValue());
    a.calcGradient(data, b.getValue() * pow(a.getValue(), b.getValue() - 1.0));
    if (a.getValue() > 0.0) {
      b.calcGradient(data, log(a.getValue()) * result);
    } else {
      b.calcGradient(data, 0.0);
    }
  }
  template<typename Data, typename Real, typename A, typename B> inline void derv11M_Pow(Data& data, const A& a, const B& b, const Real& result, const Real& multiplier) {
    checkArgumentsPow(a.getValue());
    a.calcGradient(data, multiplier * b.getValue() * pow(a.getValue(), b.getValue() - 1.0));
    if (a.getValue() > 0.0) {
      b.calcGradient(data, multiplier * log(a.getValue()) * result);
    } else {
      b.calcGradient(data, 0.0);
    }
  }
  template<typename Data, typename Real, typename A> inline void derv10_Pow(Data& data, const A& a, const typename TypeTraits<Real>::PassiveReal& b, const Real& result) {
    CODI_UNUSED(result);
    a.calcGradient(data, b * pow(a.getValue(), b - 1.0));
  }
  template<typename Data, typename Real, typename A> inline void derv10M_Pow(Data& data, const A& a, const typename TypeTraits<Real>::PassiveReal& b, const Real& result, const Real& multiplier) {
    CODI_UNUSED(result);
    a.calcGradient(data, multiplier * b * pow(a.getValue(), b - 1.0));
  }
  template<typename Data, typename Real, typename B> inline void derv01_Pow(Data& data, const typename TypeTraits<Real>::PassiveReal& a, const B& b, const Real& result) {
    checkArgumentsPow(a);
    if (a > 0.0) {
      b.calcGradient(data, log(a) * result);
    } else {
      b.calcGradient(data, 0.0);
    }
  }
  template<typename Data, typename Real, typename B> inline void derv01M_Pow(Data& data, const typename TypeTraits<Real>::PassiveReal& a, const B& b, const Real& result, const Real& multiplier) {
    checkArgumentsPow(a);
    if (a > 0.0) {
      b.calcGradient(data, multiplier * log(a) * result);
    } else {
      b.calcGradient(data, 0.0);
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
  template<typename Data, typename Real, typename A, typename B> inline void derv11_Min(Data& data, const A& a, const B& b, const Real& result) {
    if(a.getValue() < b.getValue()) {
      a.calcGradient(data);
    } else {
      b.calcGradient(data);
    }
  }
  template<typename Data, typename Real, typename A, typename B> inline void derv11M_Min(Data& data, const A& a, const B& b, const Real& result, const Real& multiplier) {
    if(a.getValue() < b.getValue()) {
      a.calcGradient(data, multiplier);
    } else {
      b.calcGradient(data, multiplier);
    }
  }
  template<typename Data, typename Real, typename A> inline void derv10_Min(Data& data, const A& a, const typename TypeTraits<Real>::PassiveReal& b, const Real& result) {
    if(a.getValue() < b) {
      a.calcGradient(data);
    }
  }
  template<typename Data, typename Real, typename A> inline void derv10M_Min(Data& data, const A& a, const typename TypeTraits<Real>::PassiveReal& b, const Real& result, const Real& multiplier) {
    if(a.getValue() < b) {
      a.calcGradient(data, multiplier);
    }
  }
  template<typename Data, typename Real, typename B> inline void derv01_Min(Data& data, const typename TypeTraits<Real>::PassiveReal& a, const B& b, const Real& result) {
    if(a >= b.getValue()) {
      b.calcGradient(data);
    }
  }
  template<typename Data, typename Real, typename B> inline void derv01M_Min(Data& data, const typename TypeTraits<Real>::PassiveReal& a, const B& b, const Real& result, const Real& multiplier) {
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
   * Implementation for f(a,b) = Max(a,b)
   */
  template<typename Data, typename Real, typename A, typename B> inline void derv11_Max(Data& data, const A& a, const B& b, const Real& result) {
    if(a.getValue() > b.getValue()) {
      a.calcGradient(data);
    } else {
      b.calcGradient(data);
    }
  }
  template<typename Data, typename Real, typename A, typename B> inline void derv11M_Max(Data& data, const A& a, const B& b, const Real& result, const Real& multiplier) {
    if(a.getValue() > b.getValue()) {
      a.calcGradient(data, multiplier);
    } else {
      b.calcGradient(data, multiplier);
    }
  }
  template<typename Data, typename Real, typename A> inline void derv10_Max(Data& data, const A& a, const typename TypeTraits<Real>::PassiveReal& b, const Real& result) {
    if(a.getValue() > b) {
      a.calcGradient(data);
    }
  }
  template<typename Data, typename Real, typename A> inline void derv10M_Max(Data& data, const A& a, const typename TypeTraits<Real>::PassiveReal& b, const Real& result, const Real& multiplier) {
    if(a.getValue() > b) {
      a.calcGradient(data, multiplier);
    }
  }
  template<typename Data, typename Real, typename B> inline void derv01_Max(Data& data, const typename TypeTraits<Real>::PassiveReal& a, const B& b, const Real& result) {
    if(a <= b.getValue()) {
      b.calcGradient(data);
    }
  }
  template<typename Data, typename Real, typename B> inline void derv01M_Max(Data& data, const typename TypeTraits<Real>::PassiveReal& a, const B& b, const Real& result, const Real& multiplier) {
    if(a <= b.getValue()) {
      b.calcGradient(data, multiplier);
    }
  }
  using std::max;
  #define NAME Max
  #define FUNCTION max
  #define PRIMAL_FUNCTION max
  #include "binaryExpression.tpp"

  #undef CODI_OPERATOR_HELPER

  /*
   * Conditional operators should behave exactly the same as with
   * non-active arguments so in each of the cases below the getValue()
   * function is called to extract the value of the expression
   */
  #define CODI_DEFINE_CONDITIONAL(OPERATOR, OP) \
    /** @brief Overload for OP with the CoDiPack expressions. @param[in] a The first argument of the operation. @param[in] b The second argument of the operation. @return The operation returns the same value the same version with double arguments. @tparam Real The real type used in the active types. @tparam A The expression for the first argument of the function @tparam B The expression for the second argument of the function*/ \
    template<typename Real, class A, class B> \
    inline bool OPERATOR(const Expression<Real, A>& a, const Expression<Real, B>& b) { \
      return a.getValue() OP b.getValue(); \
    } \
    \
    /** @brief Overload for OP with the CoDiPack expressions. @param[in] a The first argument of the operation. @param[in] b The second argument of the operation. @return The operation returns the same value the same version with double arguments. @tparam Real The real type used in the active types. @tparam A The expression for the first argument of the function */ \
    template<typename Real, class A> \
    inline bool OPERATOR(const Expression<Real, A>& a, const typename TypeTraits<Real>::PassiveReal& b) { \
      return a.getValue() OP b; \
    } \
    \
    /** @brief Overload for OP with the CoDiPack expressions. @param[in] a The first argument of the operation. @param[in] b The second argument of the operation. @return The operation returns the same value the same version with double arguments. @tparam Real The real type used in the active types. @tparam B The expression for the second argument of the function*/ \
    template<typename Real, class B> \
    inline bool OPERATOR(const typename TypeTraits<Real>::PassiveReal& a, const Expression<Real, B>& b) { \
      return a OP b.getValue(); \
    } \
    /** @brief Overload for OP with the CoDiPack expressions. @param[in] a The first argument of the operation. @param[in] b The second argument of the operation. @return The operation returns the same value the same version with double arguments. @tparam Real The real type used in the active types. @tparam A The expression for the first argument of the function */ \
    template<typename Real, class A> \
    inline bool OPERATOR(const Expression<Real, A>& a, const int& b) { \
      return a.getValue() OP b; \
    } \
    \
    /** @brief Overload for OP with the CoDiPack expressions. @param[in] a The first argument of the operation. @param[in] b The second argument of the operation. @return The operation returns the same value the same version with double arguments. @tparam Real The real type used in the active types. @tparam B The expression for the second argument of the function*/ \
    template<typename Real, class B>            \
    inline bool OPERATOR(const int& a, const Expression<Real, B>& b) { \
      return a OP b.getValue(); \
    }

  CODI_DEFINE_CONDITIONAL(operator==, ==)
  CODI_DEFINE_CONDITIONAL(operator!=, !=)
  CODI_DEFINE_CONDITIONAL(operator>, >)
  CODI_DEFINE_CONDITIONAL(operator<, <)
  CODI_DEFINE_CONDITIONAL(operator>=, >=)
  CODI_DEFINE_CONDITIONAL(operator<=, <=)

  #undef CODI_DEFINE_CONDITIONAL

  #define CODI_OPERATOR_HELPER(NAME, OP) \
    /** @brief Helper function to call operators as a function. @param[in] a The argument of the operation. @return The value of OP b @tparam A The expression for the argument of the function*/ \
    template<typename A> \
    inline auto primal_ ## NAME(const A& a) -> decltype(OP a) { \
      return OP a; \
    }

  /*
   * Now all unary functions are implemented.
   *
   * We use the naming scheme gradName for the gradient computation of the function.
   *
   */

  template<typename Real> inline typename TypeTraits<Real>::PassiveReal gradUnaryMinus(const Real& a, const Real& result) {
    CODI_UNUSED(a);
    CODI_UNUSED(result);
    return -1.0;
  }
  CODI_OPERATOR_HELPER(UnaryMinus, -)
  #define NAME UnaryMinus
  #define FUNCTION operator -
  #define PRIMAL_FUNCTION primal_UnaryMinus
  #include "unaryExpression.tpp"

  template<typename Real> inline Real gradSqrt(const Real& a, const Real& result) {
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

  template<typename Real> inline Real gradTanh(const Real& a, const Real& result) {
    CODI_UNUSED(result);
    return 1 - result * result;
  }
  using std::tanh;
  #define NAME Tanh
  #define FUNCTION tanh
  #define PRIMAL_FUNCTION tanh
  #include "unaryExpression.tpp"

  template<typename Real> inline Real gradLog(const Real& a, const Real& result) {
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

  template<typename Real> inline Real gradLog10(const Real& a, const Real& result) {
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

  template<typename Real> inline Real gradSin(const Real& a, const Real& result) {
    CODI_UNUSED(result);
    return cos(a);
  }
  using std::sin;
  #define NAME Sin
  #define FUNCTION sin
  #define PRIMAL_FUNCTION sin
  #include "unaryExpression.tpp"

  template<typename Real> inline Real gradCos(const Real& a, const Real& result) {
    CODI_UNUSED(result);
    return -sin(a);
  }
  using std::cos;
  #define NAME Cos
  #define FUNCTION cos
  #define PRIMAL_FUNCTION cos
  #include "unaryExpression.tpp"

  using std::sqrt;
  template<typename Real> inline Real gradAsin(const Real& a, const Real& result) {
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

  template<typename Real> inline Real gradAcos(const Real& a, const Real& result) {
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

  template<typename Real> inline Real gradAtan(const Real& a, const Real& result) {
    CODI_UNUSED(result);
    return 1.0 / (1 + a * a);
  }
  using std::atan;
  #define NAME Atan
  #define FUNCTION atan
  #define PRIMAL_FUNCTION atan
  #include "unaryExpression.tpp"

  template<typename Real> inline Real gradSinh(const Real& a, const Real& result) {
    CODI_UNUSED(result);
    return cosh(a);
  }
  using std::sinh;
  #define NAME Sinh
  #define FUNCTION sinh
  #define PRIMAL_FUNCTION sinh
  #include "unaryExpression.tpp"

  template<typename Real> inline Real gradCosh(const Real& a, const Real& result) {
    CODI_UNUSED(result);
    return sinh(a);
  }
  using std::cosh;
  #define NAME Cosh
  #define FUNCTION cosh
  #define PRIMAL_FUNCTION cosh
  #include "unaryExpression.tpp"

  template<typename Real> inline Real gradExp(const Real& a, const Real& result) {
    CODI_UNUSED(a);
    return result;
  }
  using std::exp;
  #define NAME Exp
  #define FUNCTION exp
  #define PRIMAL_FUNCTION exp
  #include "unaryExpression.tpp"

  template<typename Real> inline Real gradAtanh(const Real& a, const Real& result) {
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

  template<typename Real> inline Real gradAbs(const Real& a, const Real& result) {
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

  template<typename Real> inline Real gradTan(const Real& a, const Real& result) {
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
  inline Abs<Real, A> fabs(const codi::Expression<Real, A>& a) {
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
  inline const Expression<Real, A>& operator+(const Expression<Real, A>& a) {
    return a;
  }

  /***************************************************************************************
   * Functions that do not need derivatives.
   ****************************************************************************************/
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
  inline bool isinf(const codi::Expression<Real, A>& a) {
    return isinf(a.getValue());
  }

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
  inline bool isnan(const codi::Expression<Real, A>& a) {
    return isnan(a.getValue());
  }

  using std::isfinite;
  /**
   * @brief Overload for the isfinite function with expressions.
   *
   * @param[in] a The argument of the function.
   *
   * @return The result of isfinite on the primal value.
   *
   * @tparam Real The real type used in the active types.
   * @tparam A The expression for the argument of the function
   */
  template<typename Real, class A>
  inline bool isfinite(const codi::Expression<Real, A>& a) {
    return isfinite(a.getValue());
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
  inline typename codi::TypeTraits<Real>::PassiveReal floor(const codi::Expression<Real, A>& a) {
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
  inline typename codi::TypeTraits<Real>::PassiveReal ceil(const codi::Expression<Real, A>& a) {
    return ceil(a.getValue());
  }

}
