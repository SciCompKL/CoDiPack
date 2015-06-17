/**
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
#include "typeTraits.hpp"

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
     * @param[inout] gradient A helper value for forward implementations. The value is the gradient of the
     *                        lhs of the expression.
     */
    inline void calcGradient(Real& gradient) const {
      cast().calcGradient(gradient);
    }

    /**
     * @brief Calculate the gradient of the expression.
     *
     * Calculate the gradient of the mathematical operation that this
     * expression represents and pass the result to its argument.
     * For functions f(a), pass multiplier * df/da to the argument.
     *
     * @param[inout] gradient A helper value for forward implementations. The value is the gradient of the
     *                        lhs of the expression.
     * @param[in]  multiplier The Jacobi from the expression where this expression was used as an argument.
     */
    inline void calcGradient(Real& gradient, const Real& multiplier) const {
      cast().calcGradient(gradient, multiplier);
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

  /**
   * Now define particular types of expression, using static
   * polymorphism via the Curiously Recurring Template Pattern
   */

  /*
   * Enable mathematical functions with two arguments.
   *
   * @param                  OP   The name of the class
   * @param                FUNC   The name of the functions. Can also be an operator.
   * @param         PRIMAL_CALL   A function which calls the function or operator.
   * @param  DERIVATIVE_FUNC_11   The function which computes the derivative when both variables are active.
   * @param DERIVATIVE_FUNC_11M   The function which computes the derivative when both variables are active and a multiplier is given.
   * @param  DERIVATIVE_FUNC_10   The function which computes the derivative when the first variable is active.
   * @param DERIVATIVE_FUNC_10M   The function which computes the derivative when the first variable is active and a multiplier is given.
   * @param  DERIVATIVE_FUNC_01   The function which computes the derivative when the second variable is active.
   * @param DERIVATIVE_FUNC_01M   The function which computes the derivative when the second variable is active and a multiplier is given.
   */
  #define CODI_DEFINE_BINARY_FUNCTION(OP, FUNC, PRIMAL_CALL, DERIVATIVE_FUNC_11, DERIVATIVE_FUNC_11M, DERIVATIVE_FUNC_10, DERIVATIVE_FUNC_10M, DERIVATIVE_FUNC_01, DERIVATIVE_FUNC_01M)  \
    /* predefine the structs and the functions for higher order derivatives */\
    template <typename Real, class A, class B> struct OP ## 11;\
    template <typename Real, class A> struct OP ## 10;\
    template <typename Real, class B> struct OP ## 01;\
    template <typename Real, class A, class B> \
    inline  OP ## 11<Real, A, B> FUNC(const codi::Expression<Real, A>& a, const codi::Expression<Real, B>& b); \
    template <typename Real, class A> \
    inline  OP ## 10<Real, A> FUNC(const codi::Expression<Real, A>& a, const typename TypeTraits<Real>::PassiveReal& b); \
    template <typename Real, class B> \
    inline  OP ## 01<Real, B> FUNC(const typename TypeTraits<Real>::PassiveReal& a, const codi::Expression<Real, B>& b); \
    \
    /** @brief Expression implementation for OP with two active variables. @tparam Real The real type used in the active types. @tparam A The expression for the first argument of the function @tparam B The expression for the second argument of the function*/ \
    template<typename Real, class A, class B> \
    struct OP ## 11: public Expression<Real, OP ## 11<Real, A, B> > { \
      private: \
        /** @brief The first argument of the function. */ \
        const A& a_; \
        /** @brief The second argument of the function. */ \
        const B& b_; \
        /** @brief The result of the function. It is always precomputed. */ \
        Real result_; \
      public: \
        /** @brief Stores both arguments and precomputes the result of the expression. @param[in] a First argument of the expression. @param[in] b Second argument of the expression.*/ \
        OP ## 11(const Expression<Real, A>& a, const Expression<Real, B>& b) : \
          a_(a.cast()), b_(b.cast()), \
          result_(PRIMAL_CALL(a.getValue(), b.getValue())) {} \
        \
        /** @brief Calculates the jacobies of the expression and hands them down to the arguments. @details For f(x,y) it calculates df/dx and df/dy and passes these values as the multipliers to the arguments. @param[inout] gradient A helper value for forward implementations. The value is the gradient of the lhs of the expression. */ \
        inline void calcGradient(Real& gradient) const { \
          DERIVATIVE_FUNC_11(gradient, a_, b_, result_); \
        } \
        \
        /** @brief Calculates the jacobies of the expression and hands them down to the arguments. @details For f(x,y) it calculates multiplier * df/dx and multiplier * df/dy and passes these values as the multipliers to the arguments. @param[inout] gradient A helper value for forward implementations. The value is the gradient of the lhs of the expression. * @param[in]  multiplier The Jacobi from the expression where this expression was used as an argument. */ \
        inline void calcGradient(Real& gradient, const Real& multiplier) const { \
          DERIVATIVE_FUNC_11M(gradient, a_, b_, result_, multiplier); \
        } \
        \
        /** @brief Return the numerical value of the expression. @return The value of the expression. */ \
        inline const Real& getValue() const { \
          return result_; \
        } \
    }; \
    \
    /** @brief Expression implementation for OP with one active variables. @tparam Real The real type used in the active types. @tparam A The expression for the first argument of the function*/ \
    template<typename Real, class A> \
    struct OP ## 10: public Expression<Real, OP ## 10<Real, A> > { \
      private: \
        typedef typename TypeTraits<Real>::PassiveReal PassiveReal; \
        const A& a_; \
        const PassiveReal& b_; \
        Real result_; \
      public: \
        /** @brief Stores both arguments and precomputes the result of the expression. @param[in] a First argument of the expression. @param[in] b Second argument of the expression.*/ \
        OP ## 10(const Expression<Real, A>& a, const PassiveReal& b) : \
          a_(a.cast()), b_(b), \
          result_(PRIMAL_CALL(a.getValue(), b)) {} \
        \
        /** @brief Calculates the jacobies of the expression and hands them down to the arguments. @details For f(x,y) it calculates df/dx passes this value as the multiplier to the argument. @param[inout] gradient A helper value for forward implementations. The value is the gradient of the lhs of the expression. */ \
        inline void calcGradient(Real& gradient) const { \
          DERIVATIVE_FUNC_10(gradient, a_, b_, result_); \
        } \
        \
        /** @brief Calculates the jacobies of the expression and hands them down to the arguments. @details For f(x,y) it calculates multiplier * df/dx and passes this value as the multiplier to the argument. @param[inout] gradient A helper value for forward implementations. The value is the gradient of the lhs of the expression. * @param[in]  multiplier The Jacobi from the expression where this expression was used as an argument. */ \
        inline void calcGradient(Real& gradient, const Real& multiplier) const { \
          DERIVATIVE_FUNC_10M(gradient, a_, b_, result_, multiplier); \
        } \
        \
        /** @brief Return the numerical value of the expression. @return The value of the expression. */ \
        inline const Real& getValue() const { \
          return result_; \
        } \
    }; \
    \
    /** @brief Expression implementation for OP with one active variables. @tparam Real The real type used in the active types. @tparam B The expression for the second argument of the function*/ \
    template<typename Real, class B> \
    struct OP ## 01 : public Expression<Real, OP ## 01<Real, B> > { \
      private: \
        typedef typename TypeTraits<Real>::PassiveReal PassiveReal; \
        const PassiveReal& a_; \
        const B& b_; \
        Real result_; \
      public: \
        /** @brief Stores both arguments and precomputes the result of the expression. @param[in] a First argument of the expression. @param[in] b Second argument of the expression.*/ \
        OP ## 01(const PassiveReal& a, const Expression<Real, B>& b) : \
          a_(a), b_(b.cast()), \
          result_(PRIMAL_CALL(a, b.getValue())) {} \
        \
        /** @brief Calculates the jacobies of the expression and hands them down to the arguments. @details For f(x,y) it calculates df/dx passes this value as the multiplier to the argument. @param[inout] gradient A helper value for forward implementations. The value is the gradient of the lhs of the expression. */ \
        inline void calcGradient(Real& gradient) const { \
          DERIVATIVE_FUNC_01(gradient, a_, b_, result_); \
        } \
        \
        /** @brief Calculates the jacobies of the expression and hands them down to the arguments. @details For f(x,y) it calculates multiplier * df/dx and passes this value as the multiplier to the argument. @param[inout] gradient A helper value for forward implementations. The value is the gradient of the lhs of the expression. * @param[in]  multiplier The Jacobi from the expression where this expression was used as an argument. */ \
        inline void calcGradient(Real& gradient, const Real& multiplier) const { \
          DERIVATIVE_FUNC_01M(gradient, a_, b_, result_, multiplier); \
        } \
        \
        /** @brief Return the numerical value of the expression. @return The value of the expression. */ \
        inline const Real& getValue() const { \
          return result_; \
        } \
    }; \
    \
    /** @brief Overload for FUNC with the CoDiPack expressions. @param[in] a The first argument of the operation. @param[in] b The second argument of the operation. @return The implementing expression OP. @tparam Real The real type used in the active types. @tparam A The expression for the first argument of the function @tparam B The expression for the second argument of the function*/ \
    template <typename Real, class A, class B> \
    inline OP ## 11<Real, A, B> FUNC(const codi::Expression<Real, A>& a, const codi::Expression<Real, B>& b) { \
      return OP ## 11<Real, A, B>(a.cast(), b.cast()); \
    } \
    /** @brief Overload for FUNC with the CoDiPack expressions. @param[in] a The first argument of the operation. @param[in] b The second argument of the operation. @return The implementing expression OP. @tparam Real The real type used in the active types. @tparam A The expression for the first argument of the function*/ \
    template <typename Real, class A> \
    inline OP ## 10<Real, A> FUNC(const codi::Expression<Real, A>& a, const typename TypeTraits<Real>::PassiveReal& b) { \
      return OP ## 10<Real, A>(a.cast(), b); \
    } \
    /** @brief Overload for FUNC with the CoDiPack expressions. @param[in] a The first argument of the operation. @param[in] b The second argument of the operation. @return The implementing expression OP. @tparam Real The real type used in the active types. @tparam B The expression for the second argument of the function*/ \
    template <typename Real, class B> \
    inline OP ## 01<Real, B> FUNC(const typename TypeTraits<Real>::PassiveReal& a, const codi::Expression<Real, B>& b) { \
      return OP ## 01<Real, B>(a, b.cast()); \
    }

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
  template<typename Real, typename A, typename B> inline void derv11_Add(Real& gradient, const A& a, const B& b, const Real& result) {
    CODI_UNUSED(result);
    a.calcGradient(gradient);
    b.calcGradient(gradient);
  }
  template<typename Real, typename A, typename B> inline void derv11M_Add(Real& gradient, const A& a, const B& b, const Real& result, const Real& multiplier) {
    a.calcGradient(gradient, multiplier);
    b.calcGradient(gradient, multiplier);
  }
  template<typename Real, typename A> inline void derv10_Add(Real& gradient, const A& a, const typename TypeTraits<Real>::PassiveReal& b, const Real& result) {
    CODI_UNUSED(result);
    a.calcGradient(gradient);
  }
  template<typename Real, typename A> inline void derv10M_Add(Real& gradient, const A& a, const typename TypeTraits<Real>::PassiveReal& b, const Real& result, const Real& multiplier) {
    CODI_UNUSED(result);
    a.calcGradient(gradient, multiplier);
  }
  template<typename Real, typename B> inline void derv01_Add(Real& gradient, const typename TypeTraits<Real>::PassiveReal& a, const B& b, const Real& result) {
    CODI_UNUSED(result);
    b.calcGradient(gradient);
  }
  template<typename Real, typename B> inline void derv01M_Add(Real& gradient, const typename TypeTraits<Real>::PassiveReal& a, const B& b, const Real& result, const Real& multiplier) {
    CODI_UNUSED(result);
    b.calcGradient(gradient, multiplier);
  }
  CODI_OPERATOR_HELPER(Add, +)
  CODI_DEFINE_BINARY_FUNCTION(Add, operator +, primal_Add, derv11_Add, derv11M_Add, derv10_Add, derv10M_Add, derv01_Add, derv01M_Add)

  /*
   * Implementation for f(a,b) = a - b
   * df/da = 1 so simply call a
   */
  template<typename Real, typename A, typename B> inline void derv11_Subtract(Real& gradient, const A& a, const B& b, const Real& result) {
    CODI_UNUSED(result);
    a.calcGradient(gradient);
    b.calcGradient(gradient, -1.0);
  }
  template<typename Real, typename A, typename B> inline void derv11M_Subtract(Real& gradient, const A& a, const B& b, const Real& result, const Real& multiplier) {
    CODI_UNUSED(result);
    a.calcGradient(gradient, multiplier);
    b.calcGradient(gradient, -multiplier);
  }
  template<typename Real, typename A> inline void derv10_Subtract(Real& gradient, const A& a, const typename TypeTraits<Real>::PassiveReal& b, const Real& result) {
    CODI_UNUSED(result);
    a.calcGradient(gradient);
  }
  template<typename Real, typename A> inline void derv10M_Subtract(Real& gradient, const A& a, const typename TypeTraits<Real>::PassiveReal& b, const Real& result, const Real& multiplier) {
    CODI_UNUSED(result);
    a.calcGradient(gradient, multiplier);
  }
  template<typename Real, typename B> inline void derv01_Subtract(Real& gradient, const typename TypeTraits<Real>::PassiveReal& a, const B& b, const Real& result) {
    CODI_UNUSED(result);
    b.calcGradient(gradient, -1.0);
  }
  template<typename Real, typename B> inline void derv01M_Subtract(Real& gradient, const typename TypeTraits<Real>::PassiveReal& a, const B& b, const Real& result, const Real& multiplier) {
    CODI_UNUSED(result);
    b.calcGradient(gradient, -multiplier);
  }
  CODI_OPERATOR_HELPER(Subtract, -)
  CODI_DEFINE_BINARY_FUNCTION(Subtract, operator -, primal_Subtract, derv11_Subtract, derv11M_Subtract, derv10_Subtract, derv10M_Subtract, derv01_Subtract, derv01M_Subtract)

  /*
   * Implementation for f(a,b) = a * b
   */
  template<typename Real, typename A, typename B> inline void derv11_Multiply(Real& gradient, const A& a, const B& b, const Real& result) {
    CODI_UNUSED(result);
    a.calcGradient(gradient, b.getValue());
    b.calcGradient(gradient, a.getValue());
  }
  template<typename Real, typename A, typename B> inline void derv11M_Multiply(Real& gradient, const A& a, const B& b, const Real& result, const Real& multiplier) {
    CODI_UNUSED(result);
    a.calcGradient(gradient, b.getValue() * multiplier);
    b.calcGradient(gradient, a.getValue() * multiplier);
  }
  template<typename Real, typename A> inline void derv10_Multiply(Real& gradient, const A& a, const typename TypeTraits<Real>::PassiveReal& b, const Real& result) {
    CODI_UNUSED(result);
    a.calcGradient(gradient, b);
  }
  template<typename Real, typename A> inline void derv10M_Multiply(Real& gradient, const A& a, const typename TypeTraits<Real>::PassiveReal& b, const Real& result, const Real& multiplier) {
    CODI_UNUSED(result);
    a.calcGradient(gradient, b * multiplier);
  }
  template<typename Real, typename B> inline void derv01_Multiply(Real& gradient, const typename TypeTraits<Real>::PassiveReal& a, const B& b, const Real& result) {
    CODI_UNUSED(result);
    b.calcGradient(gradient, a);
  }
  template<typename Real, typename B> inline void derv01M_Multiply(Real& gradient, const typename TypeTraits<Real>::PassiveReal& a, const B& b, const Real& result, const Real& multiplier) {
    CODI_UNUSED(result);
    b.calcGradient(gradient, a * multiplier);
  }
  CODI_OPERATOR_HELPER(Multiply, *)
  CODI_DEFINE_BINARY_FUNCTION(Multiply, operator *, primal_Multiply, derv11_Multiply, derv11M_Multiply, derv10_Multiply, derv10M_Multiply, derv01_Multiply, derv01M_Multiply)

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
        CODI_EXCEPTION("Devision called with devisor of zero.");
      }
    }
  }
  template<typename Real, typename A, typename B> inline void derv11_Divide(Real& gradient, const A& a, const B& b, const Real& result) {
    checkArgumentsDivide(b.getValue());
    Real one_over_b = 1.0 / b.getValue();
    a.calcGradient(gradient, one_over_b);
    b.calcGradient(gradient, -result * one_over_b);
  }
  template<typename Real, typename A, typename B> inline void derv11M_Divide(Real& gradient, const A& a, const B& b, const Real& result, const Real& multiplier) {
    checkArgumentsDivide(b.getValue());
    Real one_over_b = multiplier / b.getValue();
    a.calcGradient(gradient, one_over_b);
    b.calcGradient(gradient, -result * one_over_b);
  }
  template<typename Real, typename A> inline void derv10_Divide(Real& gradient, const A& a, const typename TypeTraits<Real>::PassiveReal& b, const Real& result) {
    checkArgumentsDivide(b);
    Real one_over_b = 1.0 / b;
    a.calcGradient(gradient, one_over_b);
  }
  template<typename Real, typename A> inline void derv10M_Divide(Real& gradient, const A& a, const typename TypeTraits<Real>::PassiveReal& b, const Real& result, const Real& multiplier) {
    checkArgumentsDivide(b);
    Real one_over_b = multiplier / b;
    a.calcGradient(gradient, one_over_b);
  }
  template<typename Real, typename B> inline void derv01_Divide(Real& gradient, const typename TypeTraits<Real>::PassiveReal& a, const B& b, const Real& result) {
    checkArgumentsDivide(b.getValue());
    Real one_over_b = 1.0 / b.getValue();
    b.calcGradient(gradient, -result * one_over_b);
  }
  template<typename Real, typename B> inline void derv01M_Divide(Real& gradient, const typename TypeTraits<Real>::PassiveReal& a, const B& b, const Real& result, const Real& multiplier) {
    checkArgumentsDivide(b.getValue());
    Real one_over_b = multiplier / b.getValue();
    b.calcGradient(gradient, -result * one_over_b);
  }
  CODI_OPERATOR_HELPER(Divide, /)
  CODI_DEFINE_BINARY_FUNCTION(Divide, operator /, primal_Divide, derv11_Divide, derv11M_Divide, derv10_Divide, derv10M_Divide, derv01_Divide, derv01M_Divide)

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
  template<typename Real, typename A, typename B> inline void derv11_Atan2(Real& gradient, const A& a, const B& b, const Real& result) {
    CODI_UNUSED(result);
    checkArgumentsAtan2(a.getValue(), b.getValue());
    Real divisor = a.getValue() * a.getValue() + b.getValue() * b.getValue();
    divisor = 1.0 / divisor;
    a.calcGradient(gradient, b.getValue() * divisor);
    b.calcGradient(gradient, -a.getValue() * divisor);
  }
  template<typename Real, typename A, typename B> inline void derv11M_Atan2(Real& gradient, const A& a, const B& b, const Real& result, const Real& multiplier) {
    CODI_UNUSED(result);
    checkArgumentsAtan2(a.getValue(), b.getValue());
    Real divisor = a.getValue() * a.getValue() + b.getValue() * b.getValue();
    divisor = 1.0 / divisor;
    a.calcGradient(gradient, multiplier * b.getValue() * divisor);
    b.calcGradient(gradient, multiplier * -a.getValue() * divisor);
  }
  template<typename Real, typename A> inline void derv10_Atan2(Real& gradient, const A& a, const typename TypeTraits<Real>::PassiveReal& b, const Real& result) {
    CODI_UNUSED(result);
    checkArgumentsAtan2(a.getValue(), b);
    Real divisor = a.getValue() * a.getValue() + b * b;
    divisor = 1.0 / divisor;
    a.calcGradient(gradient, b * divisor);
  }
  template<typename Real, typename A> inline void derv10M_Atan2(Real& gradient, const A& a, const typename TypeTraits<Real>::PassiveReal& b, const Real& result, const Real& multiplier) {
    CODI_UNUSED(result);
    checkArgumentsAtan2(a.getValue(), b);
    Real divisor = a.getValue() * a.getValue() + b * b;
    divisor = 1.0 / divisor;
    a.calcGradient(gradient, multiplier * b * divisor);
  }
  template<typename Real, typename B> inline void derv01_Atan2(Real& gradient, const typename TypeTraits<Real>::PassiveReal& a, const B& b, const Real& result) {
    CODI_UNUSED(result);
    checkArgumentsAtan2(a, b.getValue());
    Real divisor = a * a + b.getValue() * b.getValue();
    divisor = 1.0 / divisor;
    b.calcGradient(gradient, -a * divisor);
  }
  template<typename Real, typename B> inline void derv01M_Atan2(Real& gradient, const typename TypeTraits<Real>::PassiveReal& a, const B& b, const Real& result, const Real& multiplier) {
    CODI_UNUSED(result);
    checkArgumentsAtan2(a, b.getValue());
    Real divisor = a * a + b.getValue() * b.getValue();
    divisor = 1.0 / divisor;
    b.calcGradient(gradient, multiplier * -a * divisor);
  }
  using std::atan2;
  CODI_DEFINE_BINARY_FUNCTION(Atan2, atan2, atan2, derv11_Atan2, derv11M_Atan2, derv10_Atan2, derv10M_Atan2, derv01_Atan2, derv01M_Atan2)

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
        CODI_EXCEPTION("Negative base for active exponend in pow function. (Value: %0.15e)", TypeTraits<Real>::getBaseValue(a));
      }
    }
  }
  template<typename Real, typename A, typename B> inline void derv11_Pow(Real& gradient, const A& a, const B& b, const Real& result) {
    checkArgumentsPow(a.getValue());
    a.calcGradient(gradient, b.getValue() * pow(a.getValue(), b.getValue() - 1.0));
    if (a.getValue() > 0.0) {
      b.calcGradient(gradient, log(a.getValue()) * result);
    } else {
      b.calcGradient(gradient, 0.0);
    }
  }
  template<typename Real, typename A, typename B> inline void derv11M_Pow(Real& gradient, const A& a, const B& b, const Real& result, const Real& multiplier) {
    checkArgumentsPow(a.getValue());
    a.calcGradient(gradient, multiplier * b.getValue() * pow(a.getValue(), b.getValue() - 1.0));
    if (a.getValue() > 0.0) {
      b.calcGradient(gradient, multiplier * log(a.getValue()) * result);
    } else {
      b.calcGradient(gradient, 0.0);
    }
  }
  template<typename Real, typename A> inline void derv10_Pow(Real& gradient, const A& a, const typename TypeTraits<Real>::PassiveReal& b, const Real& result) {
    CODI_UNUSED(result);
    a.calcGradient(gradient, b * pow(a.getValue(), b - 1.0));
  }
  template<typename Real, typename A> inline void derv10M_Pow(Real& gradient, const A& a, const typename TypeTraits<Real>::PassiveReal& b, const Real& result, const Real& multiplier) {
    CODI_UNUSED(result);
    a.calcGradient(gradient, multiplier * b * pow(a.getValue(), b - 1.0));
  }
  template<typename Real, typename B> inline void derv01_Pow(Real& gradient, const typename TypeTraits<Real>::PassiveReal& a, const B& b, const Real& result) {
    checkArgumentsPow(a);
    if (a > 0.0) {
      b.calcGradient(gradient, log(a) * result);
    } else {
      b.calcGradient(gradient, 0.0);
    }
  }
  template<typename Real, typename B> inline void derv01M_Pow(Real& gradient, const typename TypeTraits<Real>::PassiveReal& a, const B& b, const Real& result, const Real& multiplier) {
    checkArgumentsPow(a);
    if (a > 0.0) {
      b.calcGradient(gradient, multiplier * log(a) * result);
    } else {
      b.calcGradient(gradient, 0.0);
    }
  }
  using std::pow;
  CODI_DEFINE_BINARY_FUNCTION(Pow, pow, pow, derv11_Pow, derv11M_Pow, derv10_Pow, derv10M_Pow, derv01_Pow, derv01M_Pow)

  /*
   * Implementation for f(a,b) = Min(a,b)
   */
  template<typename Real, typename A, typename B> inline void derv11_Min(Real& gradient, const A& a, const B& b, const Real& result) {
    if(a.getValue() < b.getValue()) {
      a.calcGradient(gradient);
    } else {
      b.calcGradient(gradient);
    }
  }
  template<typename Real, typename A, typename B> inline void derv11M_Min(Real& gradient, const A& a, const B& b, const Real& result, const Real& multiplier) {
    if(a.getValue() < b.getValue()) {
      a.calcGradient(gradient, multiplier);
    } else {
      b.calcGradient(gradient, multiplier);
    }
  }
  template<typename Real, typename A> inline void derv10_Min(Real& gradient, const A& a, const typename TypeTraits<Real>::PassiveReal& b, const Real& result) {
    if(a.getValue() < b) {
      a.calcGradient(gradient);
    }
  }
  template<typename Real, typename A> inline void derv10M_Min(Real& gradient, const A& a, const typename TypeTraits<Real>::PassiveReal& b, const Real& result, const Real& multiplier) {
    if(a.getValue() < b) {
      a.calcGradient(gradient, multiplier);
    }
  }
  template<typename Real, typename B> inline void derv01_Min(Real& gradient, const typename TypeTraits<Real>::PassiveReal& a, const B& b, const Real& result) {
    if(a >= b.getValue()) {
      b.calcGradient(gradient);
    }
  }
  template<typename Real, typename B> inline void derv01M_Min(Real& gradient, const typename TypeTraits<Real>::PassiveReal& a, const B& b, const Real& result, const Real& multiplier) {
    if(a >= b.getValue()) {
      b.calcGradient(gradient, multiplier);
    }
  }
  using std::min;
  CODI_DEFINE_BINARY_FUNCTION(Min, min, min, derv11_Min, derv11M_Min, derv10_Min, derv10M_Min, derv01_Min, derv01M_Min)

  /*
   * Implementation for f(a,b) = Max(a,b)
   */
  template<typename Real, typename A, typename B> inline void derv11_Max(Real& gradient, const A& a, const B& b, const Real& result) {
    if(a.getValue() > b.getValue()) {
      a.calcGradient(gradient);
    } else {
      b.calcGradient(gradient);
    }
  }
  template<typename Real, typename A, typename B> inline void derv11M_Max(Real& gradient, const A& a, const B& b, const Real& result, const Real& multiplier) {
    if(a.getValue() > b.getValue()) {
      a.calcGradient(gradient, multiplier);
    } else {
      b.calcGradient(gradient, multiplier);
    }
  }
  template<typename Real, typename A> inline void derv10_Max(Real& gradient, const A& a, const typename TypeTraits<Real>::PassiveReal& b, const Real& result) {
    if(a.getValue() > b) {
      a.calcGradient(gradient);
    }
  }
  template<typename Real, typename A> inline void derv10M_Max(Real& gradient, const A& a, const typename TypeTraits<Real>::PassiveReal& b, const Real& result, const Real& multiplier) {
    if(a.getValue() > b) {
      a.calcGradient(gradient, multiplier);
    }
  }
  template<typename Real, typename B> inline void derv01_Max(Real& gradient, const typename TypeTraits<Real>::PassiveReal& a, const B& b, const Real& result) {
    if(a <= b.getValue()) {
      b.calcGradient(gradient);
    }
  }
  template<typename Real, typename B> inline void derv01M_Max(Real& gradient, const typename TypeTraits<Real>::PassiveReal& a, const B& b, const Real& result, const Real& multiplier) {
    if(a <= b.getValue()) {
      b.calcGradient(gradient, multiplier);
    }
  }
  using std::max;
  CODI_DEFINE_BINARY_FUNCTION(Max, max, max, derv11_Max, derv11M_Max, derv10_Max, derv10M_Max, derv01_Max, derv01M_Max)

  #undef CODI_OPERATOR_HELPER
  #undef CODI_DEFINE_BINARY_FUNCTION

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

  /**
   *  @brief Implementation for unary minus operator.
   *
   * @tparam Real The real type used in the active types.
   * @tparam A    The expression for the argument of the function.
   */
  template<typename Real, class A>
  struct UnaryMinus : public Expression<Real, UnaryMinus<Real, A> > {
   private:
     const A& a_;

   public:
   /**
    * @brief The unary minus operator.
    * @param a  The argument for the operation.
    */
    UnaryMinus(const Expression<Real, A>& a)
      : a_(a.cast()) { }


   /**
    * @brief Calculate the gradient of the expression.
    *
    * f(a) = -a => df/da = -1.0
    *
    * The operation gives the jacobi -1 to the argument.
    *
    * @param[inout] gradient A helper value for forward implementations. The value is the gradient of the
    *                        lhs of the expression.
    */
    inline void calcGradient(Real& gradient) const {
      a_.calcGradient(gradient, -1.0);
    }

   /**
    * @brief Calculate the gradient of the expression.
    *
    * f(a) = -a => df/da = -1.0
    *
    * The operation gives the jacobi -multiplier to the argument.
    *
    * @param[inout] gradient A helper value for forward implementations. The value is the gradient of the
    *                        lhs of the expression.
    * @param[in] multiplier  The multiplier from the expression which used this expression as an argument.
    */
    inline void calcGradient(Real& gradient, const Real& multiplier) const {
      a_.calcGradient(gradient, -multiplier);
    }

   /**
    * @brief Return the numerical value of the expression.
    *
    * @return The value of the expression.
    */
    inline Real getValue() const {
      return -a_.getValue();
    }
  };

  /**
   * @brief Overload for operator- with the CoDiPack expressions.
   *
   * @param[in] a The argument of the expression.
   * @return The implementing expression UnaryMinus.
   *
   * @tparam Real The real type used in the active types.
   * @tparam A The expression for the argument of the function.
   */
  template<typename Real, class A>
  inline UnaryMinus<Real, A> operator-(const Expression<Real, A>& a) {
    return UnaryMinus<Real, A>(a.cast());
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

  /*
   * Enable mathematical functions with one argument.
   *
   * @param              OP   The name of the implementing class.
   * @param            FUNC   The name of the function. FUNC(a) should be a valid call.
   * @param DERIVATIVE_FUNC   The name of the function for the derivative calculation.
   */
  # define CODI_DEFINE_UNARY_FUNCTION(OP, FUNC, DERIVATIVE_FUNC)  \
    /* predefine the struct and the function for higher order derivatives */\
    template<typename Real, class A> struct OP; \
    template <typename Real, class A> \
    inline  codi:: OP<Real, A> FUNC(const codi::Expression<Real, A>& a); \
    \
    using std:: FUNC; \
    /** @brief Expression implementation for OP. @tparam Real The real type used in the active types. @tparam A The expression for the argument of the function.*/ \
    template<typename Real, class A> \
    struct OP : public Expression<Real, OP<Real, A> > { \
      private: \
        /** @brief The argument of the function. */ \
        const A& a_; \
        /** @brief The result of the function. It is always precomputed. */ \
        Real result_; \
      public: \
        /** @brief Stores the argument and precomputes the result of the expression. @param[in] a Argument of the expression.*/ \
        OP(const Expression<Real, A>& a) : \
          a_(a.cast()), \
          result_(FUNC(a.getValue())) {} \
      \
      /** @brief Calculates the jacobie of the expression and hands them down to the argument. @details For f(x) it calculates df/dx and passes this value as the multiplier to the argument. @param[inout] gradient A helper value for forward implementations. The value is the gradient of the lhs of the expression. */ \
      inline void calcGradient(Real& gradient) const { \
        a_.calcGradient(gradient, DERIVATIVE_FUNC(a_.getValue(), result_)); \
      } \
      \
      /** @brief Calculates the jacobi of the expression and hands them down to the argument. @details For f(x) it calculates multiplier * df/dx and passes this value as the multiplier to the argument. @param[inout] gradient A helper value for forward implementations. The value is the gradient of the lhs of the expression. * @param[in]  multiplier The Jacobi from the expression where this expression was used as an argument. */ \
      inline void calcGradient(Real& gradient, const Real& multiplier) const { \
        a_.calcGradient(gradient, DERIVATIVE_FUNC(a_.getValue(), result_)*multiplier); \
      } \
      \
      /** @brief Return the numerical value of the expression. @return The value of the expression. */ \
      inline const Real& getValue() const { \
        return result_; \
      } \
    }; \
    \
    /** @brief Overload for FUNC with the CoDiPack expressions. @param[in] a The argument of the operation. @return The implementing expression OP. @tparam Real The real type used in the active types. @tparam A The expression for the first argument of the function. */ \
    template <typename Real, class A> \
    inline codi:: OP<Real, A> FUNC(const codi::Expression<Real, A>& a) { \
      return codi:: OP<Real, A>(a.cast()); \
    }

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
  CODI_DEFINE_UNARY_FUNCTION(Sqrt, sqrt, gradSqrt)

  template<typename Real> inline Real gradTanh(const Real& a, const Real& result) {
    CODI_UNUSED(result);
    return 1 - result * result;
  }
  CODI_DEFINE_UNARY_FUNCTION(Tanh, tanh, gradTanh)

  template<typename Real> inline Real gradLog(const Real& a, const Real& result) {
    CODI_UNUSED(result);
    if(CheckExpressionArguments) {
      if(0.0 > TypeTraits<Real>::getBaseValue(a)) {
        CODI_EXCEPTION("Logarithm of negative value or zero.(Value: %0.15e)", TypeTraits<Real>::getBaseValue(a));
      }
    }
    return 1.0 / a;
  }
  CODI_DEFINE_UNARY_FUNCTION(Log, log, gradLog)

  template<typename Real> inline Real gradLog10(const Real& a, const Real& result) {
    CODI_UNUSED(result);
    if(CheckExpressionArguments) {
      if(0.0 > TypeTraits<Real>::getBaseValue(a)) {
        CODI_EXCEPTION("Logarithm of negative value or zero.(Value: %0.15e)", TypeTraits<Real>::getBaseValue(a));
      }
    }
    return 0.434294481903252 / a;
  }
  CODI_DEFINE_UNARY_FUNCTION(Log10, log10, gradLog10)

  template<typename Real> inline Real gradSin(const Real& a, const Real& result) {
    CODI_UNUSED(result);
    return cos(a);
  }
  CODI_DEFINE_UNARY_FUNCTION(Sin, sin, gradSin)

  template<typename Real> inline Real gradCos(const Real& a, const Real& result) {
    CODI_UNUSED(result);
    return -sin(a);
  }
  CODI_DEFINE_UNARY_FUNCTION(Cos, cos, gradCos)

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
  CODI_DEFINE_UNARY_FUNCTION(Asin, asin, gradAsin)

  template<typename Real> inline Real gradAcos(const Real& a, const Real& result) {
    CODI_UNUSED(result);
    if(CheckExpressionArguments) {
      if(TypeTraits<Real>::getBaseValue(a) <= -1.0 || 1.0 <= TypeTraits<Real>::getBaseValue(a)) {
        CODI_EXCEPTION("acos outside of (-1, 1).(Value: %0.15e)", TypeTraits<Real>::getBaseValue(a));
      }
    }
    return -1.0 / sqrt(1.0 - a * a);
  }
  CODI_DEFINE_UNARY_FUNCTION(Acos, acos, gradAcos)

  template<typename Real> inline Real gradAtan(const Real& a, const Real& result) {
    CODI_UNUSED(result);
    return 1.0 / (1 + a * a);
  }
  CODI_DEFINE_UNARY_FUNCTION(Atan, atan, gradAtan)

  template<typename Real> inline Real gradSinh(const Real& a, const Real& result) {
    CODI_UNUSED(result);
    return cosh(a);
  }
  CODI_DEFINE_UNARY_FUNCTION(Sinh, sinh, gradSinh)

  template<typename Real> inline Real gradCosh(const Real& a, const Real& result) {
    CODI_UNUSED(result);
    return sinh(a);
  }
  CODI_DEFINE_UNARY_FUNCTION(Cosh, cosh, gradCosh)

  template<typename Real> inline Real gradExp(const Real& a, const Real& result) {
    CODI_UNUSED(a);
    return result;
  }
  CODI_DEFINE_UNARY_FUNCTION(Exp, exp, gradExp)

  template<typename Real> inline Real gradAtanh(const Real& a, const Real& result) {
    CODI_UNUSED(result);
    if(CheckExpressionArguments) {
      if(TypeTraits<Real>::getBaseValue(a) <= -1.0 || 1.0 <= TypeTraits<Real>::getBaseValue(a)) {
        CODI_EXCEPTION("atanh outside of (-1, 1).(Value: %0.15e)", TypeTraits<Real>::getBaseValue(a));
      }
    }
    return 1.0 / (1 - a * a);
  }
  CODI_DEFINE_UNARY_FUNCTION(Atanh, atanh, gradAtanh)

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
  CODI_DEFINE_UNARY_FUNCTION(Abs, abs, gradAbs)

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
  CODI_DEFINE_UNARY_FUNCTION(Tan, tan, gradTan)
  # undef CODI_DEFINE_UNARY_FUNCTION


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

  /***************************************************************************************
   * Functions which do not need derivatives.
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
