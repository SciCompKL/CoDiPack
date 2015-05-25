/**
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015 Chair for Scientific Computing, TU Kaiserslautern
 *
 * This file is part of CoDiPack.
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
 * Authors: TODO
 */
#pragma once

#include <cmath>
#include <algorithm>

#include "configure.h"
#include "typeTraits.hpp"

namespace codi {

  /**
   * The Expression type from which all other types of expression
   * derive. Each member function simply calls the specialized version
   * of the function according to the expression's true type, which is
   * given by its template argument.
   *
   * @template Real  The data type of the primal values and the gradient values.
   */

  template<typename Real, class A>
  struct Expression {

    typedef typename TypeTraits<Real>::PassiveReal PassiveReal;

    /**
     * Cast the expression to its true type, given by the template
     * argument
     */
    inline const A& cast() const {
      return static_cast<const A&>(*this);
    }

    /**
     * Calculate the gradient of the mathematical operation that this
     * expression represents and pass the result to its argument.
     * For functions f(a), pass df/da to the argument in the
     * first case and pass multiplier*df/da in the second case.
     */
    inline void calcGradient(Real& gradient) const {
      cast().calcGradient(gradient);
    }

    /**
     * As the previous but multiplying the gradient by "multiplier"
     */
    inline void calcGradient(Real& gradient, const Real& multiplier) const {
      cast().calcGradient(gradient, multiplier);
    }

    /**
     * Return the numerical value of the expression
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

  /**
   * Multiply: an expression multiplied by another expression
   */

  template <typename Real, class A, class B>
  struct Multiply : public Expression<Real, Multiply<Real, A,B> > {
    Multiply(const Expression<Real, A>& a, const Expression<Real, B>& b)
      : a_(a.cast()), b_(b.cast()), result_(a_.getValue()*b_.getValue()) { }
    // If f(a,b) = a*b then df/da = b and df/db = a
    inline void calcGradient(Real& gradient) const {
      a_.calcGradient(gradient, b_.getValue());
      b_.calcGradient(gradient, a_.getValue());
    }
    inline void calcGradient(Real& gradient, const Real& multiplier) const {
      a_.calcGradient(gradient, b_.getValue()*multiplier);
      b_.calcGradient(gradient, a_.getValue()*multiplier);
    }
    inline const Real& getValue() const {
      return result_;
    }
  private:
    const A& a_;
    const B& b_;
    Real result_;
  };

  /**
   * Enable mathematical functions with two arguments.
   */
  # define CODI_DEFINE_BINARY_FUNCTION(OP, FUNC, PRIMAL_CALL, DERIVATIVE_FUNC_11, DERIVATIVE_FUNC_11M, DERIVATIVE_FUNC_10, DERIVATIVE_FUNC_10M, DERIVATIVE_FUNC_01, DERIVATIVE_FUNC_01M)  \
  /* predefine the struct and the function for higher order derivatrives */\
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
      template<typename Real, class A, class B> \
      struct OP ## 11: public Expression<Real, OP ## 11<Real, A, B> > { \
        private: \
          const A& a_; \
          const B& b_; \
          Real result_; \
        public: \
          OP ## 11(const Expression<Real, A>& a, const Expression<Real, B>& b) : \
            a_(a.cast()), b_(b.cast()), \
            result_(PRIMAL_CALL(a.getValue(), b.getValue())) {} \
        \
        inline void calcGradient(Real& gradient) const { \
          DERIVATIVE_FUNC_11(gradient, a_, b_, result_); \
        } \
        \
        inline void calcGradient(Real& gradient, const Real& multiplier) const { \
          DERIVATIVE_FUNC_11M(gradient, a_, b_, result_, multiplier); \
        } \
        \
        inline const Real& getValue() const { \
          return result_; \
        } \
      }; \
      \
      template<typename Real, class A> \
      struct OP ## 10: public Expression<Real, OP ## 10<Real, A> > { \
        public: \
          typedef typename TypeTraits<Real>::PassiveReal PassiveReal; \
        private: \
          const A& a_; \
          const PassiveReal& b_; \
          Real result_; \
        public: \
          OP ## 10(const Expression<Real, A>& a, const PassiveReal& b) : \
            a_(a.cast()), b_(b), \
            result_(PRIMAL_CALL(a.getValue(), b)) {} \
        \
        inline void calcGradient(Real& gradient) const { \
          DERIVATIVE_FUNC_10(gradient, a_, b_, result_); \
        } \
        \
        inline void calcGradient(Real& gradient, const Real& multiplier) const { \
          DERIVATIVE_FUNC_10M(gradient, a_, b_, result_, multiplier); \
        } \
        \
        inline const Real& getValue() const { \
          return result_; \
        } \
      }; \
      \
      template<typename Real, class B> \
      struct OP ## 01 : public Expression<Real, OP ## 01<Real, B> > { \
        public: \
          typedef typename TypeTraits<Real>::PassiveReal PassiveReal; \
        private: \
          const PassiveReal& a_; \
          const B& b_; \
          Real result_; \
        public: \
          OP ## 01(const PassiveReal& a, const Expression<Real, B>& b) : \
            a_(a), b_(b.cast()), \
            result_(PRIMAL_CALL(a, b.getValue())) {} \
        \
        inline void calcGradient(Real& gradient) const { \
          DERIVATIVE_FUNC_01(gradient, a_, b_, result_); \
        } \
        \
        inline void calcGradient(Real& gradient, const Real& multiplier) const { \
          DERIVATIVE_FUNC_01M(gradient, a_, b_, result_, multiplier); \
        } \
        \
        inline const Real& getValue() const { \
          return result_; \
        } \
      }; \
      \
      template <typename Real, class A, class B> \
      inline OP ## 11<Real, A, B> FUNC(const codi::Expression<Real, A>& a, const codi::Expression<Real, B>& b) { \
        return OP ## 11<Real, A, B>(a.cast(), b.cast()); \
      } \
      template <typename Real, class A> \
      inline OP ## 10<Real, A> FUNC(const codi::Expression<Real, A>& a, const typename TypeTraits<Real>::PassiveReal& b) { \
        return OP ## 10<Real, A>(a.cast(), b); \
      } \
      template <typename Real, class B> \
      inline OP ## 01<Real, B> FUNC(const typename TypeTraits<Real>::PassiveReal& a, const codi::Expression<Real, B>& b) { \
        return OP ## 01<Real, B>(a, b.cast()); \
      }

  #define CODI_OPERATOR_HELPER(NAME, OP) \
    template<typename A, typename B> \
    inline auto primal_ ## NAME(const A& a, const B& b) -> decltype(a OP b) { \
      return a OP b; \
    }

  /**
   * If f(a,b) = a + b, df/da = 1 and likewise for df/db so simply
   * call a and b's versions of calcGradient
   */
  template<typename Real, typename A, typename B> inline void derv11_Add(Real& gradient, const A& a, const B& b, const Real& /*result*/) {
    a.calcGradient(gradient);
    b.calcGradient(gradient);
  }
  template<typename Real, typename A, typename B> inline void derv11M_Add(Real& gradient, const A& a, const B& b, const Real& /*result*/, const Real& multiplier) {
    a.calcGradient(gradient, multiplier);
    b.calcGradient(gradient, multiplier);
  }
  template<typename Real, typename A> inline void derv10_Add(Real& gradient, const A& a, const typename TypeTraits<Real>::PassiveReal& b, const Real& /*result*/) {
    a.calcGradient(gradient);
  }
  template<typename Real, typename A> inline void derv10M_Add(Real& gradient, const A& a, const typename TypeTraits<Real>::PassiveReal& b, const Real& /*result*/, const Real& multiplier) {
    a.calcGradient(gradient, multiplier);
  }
  template<typename Real, typename B> inline void derv01_Add(Real& gradient, const typename TypeTraits<Real>::PassiveReal& a, const B& b, const Real& /*result*/) {
    b.calcGradient(gradient);
  }
  template<typename Real, typename B> inline void derv01M_Add(Real& gradient, const typename TypeTraits<Real>::PassiveReal& a, const B& b, const Real& /*result*/, const Real& multiplier) {
    b.calcGradient(gradient, multiplier);
  }
  CODI_OPERATOR_HELPER(Add, +)
  CODI_DEFINE_BINARY_FUNCTION(Add, operator +, primal_Add, derv11_Add, derv11M_Add, derv10_Add, derv10M_Add, derv01_Add, derv01M_Add)

  /**
   * If f(a,b) = a - b, df/da = 1 so simply
   * call a
   */
  template<typename Real, typename A, typename B> inline void derv11_Subtract(Real& gradient, const A& a, const B& b, const Real& /*result*/) {
    a.calcGradient(gradient);
    b.calcGradient(gradient, -1.0);
  }
  template<typename Real, typename A, typename B> inline void derv11M_Subtract(Real& gradient, const A& a, const B& b, const Real& /*result*/, const Real& multiplier) {
    a.calcGradient(gradient, multiplier);
    b.calcGradient(gradient, -multiplier);
  }
  template<typename Real, typename A> inline void derv10_Subtract(Real& gradient, const A& a, const typename TypeTraits<Real>::PassiveReal& b, const Real& /*result*/) {
    a.calcGradient(gradient);
  }
  template<typename Real, typename A> inline void derv10M_Subtract(Real& gradient, const A& a, const typename TypeTraits<Real>::PassiveReal& b, const Real& /*result*/, const Real& multiplier) {
    a.calcGradient(gradient, multiplier);
  }
  template<typename Real, typename B> inline void derv01_Subtract(Real& gradient, const typename TypeTraits<Real>::PassiveReal& a, const B& b, const Real& /*result*/) {
    b.calcGradient(gradient, -1.0);
  }
  template<typename Real, typename B> inline void derv01M_Subtract(Real& gradient, const typename TypeTraits<Real>::PassiveReal& a, const B& b, const Real& /*result*/, const Real& multiplier) {
    b.calcGradient(gradient, -multiplier);
  }
  CODI_OPERATOR_HELPER(Subtract, -)
  CODI_DEFINE_BINARY_FUNCTION(Subtract, operator -, primal_Subtract, derv11_Subtract, derv11M_Subtract, derv10_Subtract, derv10M_Subtract, derv01_Subtract, derv01M_Subtract)

  template<typename Real, typename A, typename B> inline void derv11_Multiply(Real& gradient, const A& a, const B& b, const Real& /*result*/) {
    a.calcGradient(gradient, b.getValue());
    b.calcGradient(gradient, a.getValue());
  }
  template<typename Real, typename A, typename B> inline void derv11M_Multiply(Real& gradient, const A& a, const B& b, const Real& /*result*/, const Real& multiplier) {
    a.calcGradient(gradient, b.getValue() * multiplier);
    b.calcGradient(gradient, a.getValue() * multiplier);
  }
  template<typename Real, typename A> inline void derv10_Multiply(Real& gradient, const A& a, const typename TypeTraits<Real>::PassiveReal& b, const Real& /*result*/) {
    a.calcGradient(gradient, b);
  }
  template<typename Real, typename A> inline void derv10M_Multiply(Real& gradient, const A& a, const typename TypeTraits<Real>::PassiveReal& b, const Real& /*result*/, const Real& multiplier) {
    a.calcGradient(gradient, b * multiplier);
  }
  template<typename Real, typename B> inline void derv01_Multiply(Real& gradient, const typename TypeTraits<Real>::PassiveReal& a, const B& b, const Real& /*result*/) {
    b.calcGradient(gradient, a);
  }
  template<typename Real, typename B> inline void derv01M_Multiply(Real& gradient, const typename TypeTraits<Real>::PassiveReal& a, const B& b, const Real& /*result*/, const Real& multiplier) {
    b.calcGradient(gradient, a * multiplier);
  }
  CODI_OPERATOR_HELPER(Multiply, *)
  CODI_DEFINE_BINARY_FUNCTION(Multiply, operator *, primal_Multiply, derv11_Multiply, derv11M_Multiply, derv10_Multiply, derv10M_Multiply, derv01_Multiply, derv01M_Multiply)

  template<typename Real, typename A, typename B> inline void derv11_Divide(Real& gradient, const A& a, const B& b, const Real& result) {
    Real one_over_b = 1.0 / b.getValue();
    a.calcGradient(gradient, one_over_b);
    b.calcGradient(gradient, -result * one_over_b);
  }
  template<typename Real, typename A, typename B> inline void derv11M_Divide(Real& gradient, const A& a, const B& b, const Real& result, const Real& multiplier) {
    Real one_over_b = multiplier / b.getValue();
    a.calcGradient(gradient, one_over_b);
    b.calcGradient(gradient, -result * one_over_b);
  }
  template<typename Real, typename A> inline void derv10_Divide(Real& gradient, const A& a, const typename TypeTraits<Real>::PassiveReal& b, const Real& result) {
    Real one_over_b = 1.0 / b;
    a.calcGradient(gradient, one_over_b);
  }
  template<typename Real, typename A> inline void derv10M_Divide(Real& gradient, const A& a, const typename TypeTraits<Real>::PassiveReal& b, const Real& result, const Real& multiplier) {
    Real one_over_b = multiplier / b;
    a.calcGradient(gradient, one_over_b);
  }
  template<typename Real, typename B> inline void derv01_Divide(Real& gradient, const typename TypeTraits<Real>::PassiveReal& a, const B& b, const Real& result) {
    Real one_over_b = 1.0 / b.getValue();
    b.calcGradient(gradient, -result * one_over_b);
  }
  template<typename Real, typename B> inline void derv01M_Divide(Real& gradient, const typename TypeTraits<Real>::PassiveReal& a, const B& b, const Real& result, const Real& multiplier) {
    Real one_over_b = multiplier / b.getValue();
    b.calcGradient(gradient, -result * one_over_b);
  }
  CODI_OPERATOR_HELPER(Divide, /)
  CODI_DEFINE_BINARY_FUNCTION(Divide, operator /, primal_Divide, derv11_Divide, derv11M_Divide, derv10_Divide, derv10M_Divide, derv01_Divide, derv01M_Divide)

  template<typename Real, typename A, typename B> inline void derv11_Atan2(Real& gradient, const A& a, const B& b, const Real& /*result*/) {
    Real divisor = a.getValue() * a.getValue() + b.getValue() * b.getValue();
    divisor = 1.0 / divisor;
    a.calcGradient(gradient, b.getValue() * divisor);
    b.calcGradient(gradient, -a.getValue() * divisor);
  }
  template<typename Real, typename A, typename B> inline void derv11M_Atan2(Real& gradient, const A& a, const B& b, const Real& /*result*/, const Real& multiplier) {
    Real divisor = a.getValue() * a.getValue() + b.getValue() * b.getValue();
    divisor = 1.0 / divisor;
    a.calcGradient(gradient, multiplier * b.getValue() * divisor);
    b.calcGradient(gradient, multiplier * -a.getValue() * divisor);
  }
  template<typename Real, typename A> inline void derv10_Atan2(Real& gradient, const A& a, const typename TypeTraits<Real>::PassiveReal& b, const Real& /*result*/) {
    Real divisor = a.getValue() * a.getValue() + b * b;
    divisor = 1.0 / divisor;
    a.calcGradient(gradient, b * divisor);
  }
  template<typename Real, typename A> inline void derv10M_Atan2(Real& gradient, const A& a, const typename TypeTraits<Real>::PassiveReal& b, const Real& /*result*/, const Real& multiplier) {
    Real divisor = a.getValue() * a.getValue() + b * b;
    divisor = 1.0 / divisor;
    a.calcGradient(gradient, multiplier * b * divisor);
  }
  template<typename Real, typename B> inline void derv01_Atan2(Real& gradient, const typename TypeTraits<Real>::PassiveReal& a, const B& b, const Real& /*result*/) {
    Real divisor = a * a + b.getValue() * b.getValue();
    divisor = 1.0 / divisor;
    b.calcGradient(gradient, -a * divisor);
  }
  template<typename Real, typename B> inline void derv01M_Atan2(Real& gradient, const typename TypeTraits<Real>::PassiveReal& a, const B& b, const Real& /*result*/, const Real& multiplier) {
    Real divisor = a * a + b.getValue() * b.getValue();
    divisor = 1.0 / divisor;
    b.calcGradient(gradient, multiplier * -a * divisor);
  }
  using std::atan2;
  CODI_DEFINE_BINARY_FUNCTION(Atan2, atan2, atan2, derv11_Atan2, derv11M_Atan2, derv10_Atan2, derv10M_Atan2, derv01_Atan2, derv01M_Atan2)

  template<typename Real, typename A, typename B> inline void derv11_Pow(Real& gradient, const A& a, const B& b, const Real& result) {
    a.calcGradient(gradient, b.getValue() * pow(a.getValue(), b.getValue() - 1.0));
    if (a.getValue() > 0.0) {
      return b.calcGradient(gradient, log(a.getValue()) * result);
    } else {
      return b.calcGradient(gradient, 0.0);
    }
  }
  template<typename Real, typename A, typename B> inline void derv11M_Pow(Real& gradient, const A& a, const B& b, const Real& result, const Real& multiplier) {
    a.calcGradient(gradient, multiplier * b.getValue() * pow(a.getValue(), b.getValue() - 1.0));
    if (a.getValue() > 0.0) {
      return b.calcGradient(gradient, multiplier * log(a.getValue()) * result);
    } else {
      return b.calcGradient(gradient, 0.0);
    }
  }
  template<typename Real, typename A> inline void derv10_Pow(Real& gradient, const A& a, const typename TypeTraits<Real>::PassiveReal& b, const Real& /*result*/) {
    a.calcGradient(gradient, b * pow(a.getValue(), b - 1.0));
  }
  template<typename Real, typename A> inline void derv10M_Pow(Real& gradient, const A& a, const typename TypeTraits<Real>::PassiveReal& b, const Real& /*result*/, const Real& multiplier) {
    a.calcGradient(gradient, multiplier * b * pow(a.getValue(), b - 1.0));
  }
  template<typename Real, typename B> inline void derv01_Pow(Real& gradient, const typename TypeTraits<Real>::PassiveReal& a, const B& b, const Real& result) {
    if (a > 0.0) {
      return b.calcGradient(gradient, log(a) * result);
    } else {
      return b.calcGradient(gradient, 0.0);
    }
  }
  template<typename Real, typename B> inline void derv01M_Pow(Real& gradient, const typename TypeTraits<Real>::PassiveReal& a, const B& b, const Real& result, const Real& multiplier) {
    if (a > 0.0) {
      return b.calcGradient(gradient, multiplier * log(a) * result);
    } else {
      return b.calcGradient(gradient, 0.0);
    }
  }
  using std::pow;
  CODI_DEFINE_BINARY_FUNCTION(Pow, pow, pow, derv11_Pow, derv11M_Pow, derv10_Pow, derv10M_Pow, derv01_Pow, derv01M_Pow)

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
#undef CODI_OPERATOR_HELPER;
#undef CODI_DEFINE_BINARY_FUNCTION;

  /**
   *  Conditional operators should behave exactly the same as with
   * non-active arguments so in each of the cases below the getValue()
   * function is called to extract the value of the expression
   */
#define CODI_DEFINE_CONDITIONAL(OPERATOR, OP)      \
template<typename Real, class A, class B>          \
inline              \
bool OPERATOR(const Expression<Real, A>& a,        \
      const Expression<Real, B>& b) {      \
  return a.getValue() OP b.getValue();        \
}                \
                              \
template<typename Real, class A>            \
inline              \
bool OPERATOR(const Expression<Real, A>& a, const Real& b) {  \
  return a.getValue() OP b;          \
}                \
                              \
template<typename Real, class B>            \
inline              \
bool OPERATOR(const Real& a, const Expression<Real, B>& b) {  \
  return a OP b.getValue();          \
}\
  template<typename Real, class A>            \
  inline              \
  bool OPERATOR(const Expression<Real, A>& a, const int& b) {  \
    return a.getValue() OP b;          \
  }                \
                                \
  template<typename Real, class B>            \
  inline              \
  bool OPERATOR(const int& a, const Expression<Real, B>& b) {  \
    return a OP b.getValue();          \
  }
  CODI_DEFINE_CONDITIONAL(operator==, ==)

  CODI_DEFINE_CONDITIONAL(operator!=, !=)

  CODI_DEFINE_CONDITIONAL(operator>, >)

  CODI_DEFINE_CONDITIONAL(operator<, <)

  CODI_DEFINE_CONDITIONAL(operator>=, >=)

  CODI_DEFINE_CONDITIONAL(operator<=, <=)

#undef CODI_DEFINE_CONDITIONAL

  /**
   *  UnaryMinus: negation of expression
   */
  template<typename Real, class A>
  struct UnaryMinus : public Expression<Real, UnaryMinus<Real, A> > {
    UnaryMinus(const Expression<Real, A>& a)
      : a_(a.cast()) { }

    inline void calcGradient(Real& gradient) const {
      a_.calcGradient(gradient, -1.0);
    }

    inline void calcGradient(Real& gradient, const Real& multiplier) const {
      a_.calcGradient(gradient, -multiplier);
    }

    inline Real getValue() const {
      return -a_.getValue();
    }

  private:
    const A& a_;
  };

  /**
   *  Overload unary minus of expression
   */
  template<typename Real, class A>
  inline
  UnaryMinus<Real, A> operator-(const Expression<Real, A>& a) {
    return UnaryMinus<Real, A>(a.cast());
  }

  /**
   *  Unary plus: returns the argument
   */
  template<typename Real, class A>
  inline
  const Expression<Real, A>& operator+(const Expression<Real, A>& a) {
    return a;
  }

/**
 * Enable mathematical functions with one argument.
 */
# define CODI_DEFINE_UNARY_FUNCTION(OP, FUNC, DERIVATIVE_FUNC)  \
/* predefine the struct and the function for higher order derivatrives */\
    template<typename Real, class A> struct OP; \
    template <typename Real, class A> \
    inline  codi:: OP<Real, A> FUNC(const codi::Expression<Real, A>& a); \
    \
    using std:: FUNC; \
    template<typename Real, class A> \
    struct OP : public Expression<Real, OP<Real, A> > { \
      private: \
        const A& a_; \
        Real result_; \
      public: \
        OP(const Expression<Real, A>& a) : \
          a_(a.cast()), \
          result_(FUNC(a.getValue())) {} \
      \
      inline void calcGradient(Real& gradient) const { \
        a_.calcGradient(gradient, DERIVATIVE_FUNC(a_.getValue(), result_)); \
      } \
      \
      inline void calcGradient(Real& gradient, const Real& multiplier) const { \
        a_.calcGradient(gradient, DERIVATIVE_FUNC(a_.getValue(), result_)*multiplier); \
      } \
      \
      inline const Real& getValue() const { \
        return result_; \
      } \
    }; \
    \
    template <typename Real, class A> \
    inline codi:: OP<Real, A> FUNC(const codi::Expression<Real, A>& a) { \
      return codi:: OP<Real, A>(a.cast()); \
    }

    template<typename Real> inline Real gradSqrt(const Real& /*a*/, const Real& result) {
      if(result != 0.0) {
        return 0.5 / result;
      } else {
        return (Real)0.0;
      }
    }
    CODI_DEFINE_UNARY_FUNCTION(Sqrt, sqrt, gradSqrt)

    template<typename Real> inline Real gradTanh(const Real& /*a*/, const Real& result) {
      return 1 - result * result;
    }
    CODI_DEFINE_UNARY_FUNCTION(Tanh, tanh, gradTanh)

    template<typename Real> inline Real gradLog(const Real& a, const Real& /*result*/) {
      return 1.0 / a;
    }
    CODI_DEFINE_UNARY_FUNCTION(Log, log, gradLog)

    template<typename Real> inline Real gradLog10(const Real& a, const Real& /*result*/) {
      return 0.434294481903252 / a;
    }
    CODI_DEFINE_UNARY_FUNCTION(Log10, log10, gradLog10)

    template<typename Real> inline Real gradSin(const Real& a, const Real& /*result*/) {
      return cos(a);
    }
    CODI_DEFINE_UNARY_FUNCTION(Sin, sin, gradSin)

    template<typename Real> inline Real gradCos(const Real& a, const Real& /*result*/) {
      return -sin(a);
    }
    CODI_DEFINE_UNARY_FUNCTION(Cos, cos, gradCos)

    using std::sqrt;
    template<typename Real> inline Real gradAsin(const Real& a, const Real& /*result*/) {
      return 1.0 / sqrt(1.0 - a * a);
    }
    CODI_DEFINE_UNARY_FUNCTION(Asin, asin, gradAsin)

    template<typename Real> inline Real gradAcos(const Real& a, const Real& /*result*/) {
      return -1.0 / sqrt(1.0 - a * a);
    }
    CODI_DEFINE_UNARY_FUNCTION(Acos, acos, gradAcos)

    template<typename Real> inline Real gradAtan(const Real& a, const Real& /*result*/) {
      return 1.0 / (1 + a * a);
    }
    CODI_DEFINE_UNARY_FUNCTION(Atan, atan, gradAtan)

    template<typename Real> inline Real gradSinh(const Real& a, const Real& /*result*/) {
      return cosh(a);
    }
    CODI_DEFINE_UNARY_FUNCTION(Sinh, sinh, gradSinh)

    template<typename Real> inline Real gradCosh(const Real& a, const Real& /*result*/) {
      return sinh(a);
    }
    CODI_DEFINE_UNARY_FUNCTION(Cosh, cosh, gradCosh)

    template<typename Real> inline Real gradExp(const Real& /*a*/, const Real& result) {
      return result;
    }
    CODI_DEFINE_UNARY_FUNCTION(Exp, exp, gradExp)

    template<typename Real> inline Real gradAtanh(const Real& a, const Real& /*result*/) {
      return 1.0 / (1 - a * a);
    }
    CODI_DEFINE_UNARY_FUNCTION(Atanh, atanh, gradAtanh)

    template<typename Real> inline Real gradAbs(const Real& a, const Real& /*result*/) {
      if(a < 0.0) {
        return (Real)-1.0;
      } else if(a > 0.0) {
        return (Real)1.0;
      } else {
        return (Real)0.0;
      }
    }
    CODI_DEFINE_UNARY_FUNCTION(Abs, abs, gradAbs)

    template<typename Real> inline Real gradTan(const Real& a, const Real& /*result*/) {
      Real tmp = 1 / cos(a);
      return tmp * tmp;
    }
    CODI_DEFINE_UNARY_FUNCTION(Tan, tan, gradTan)
# undef CODI_DEFINE_UNARY_FUNCTION

template<typename Real, class A>
inline Abs<Real, A> fabs(const codi::Expression<Real, A>& a) {
  return Abs<Real, A>(a.cast());
}

/**
 * Need to add ceil, floor...
 * Lots more math function in math.h: erf, bessel functions etc...
 */
template<typename Real, class A>
inline bool isinf(const codi::Expression<Real, A>& a) {
  return isinf(a.getValue());
}

template<typename Real, class A>
inline bool isnan(const codi::Expression<Real, A>& a) {
  return isnan(a.getValue());
}

using std::isfinite;
template<typename Real, class A>
inline bool isfinite(const codi::Expression<Real, A>& a) {
  return isfinite(a.getValue());
}

using std::floor;
template<typename Real, class A>
inline typename codi::TypeTraits<Real>::PassiveReal floor(const codi::Expression<Real, A>& a) {
  return floor(a.getValue());
}

using std::ceil;
template<typename Real, class A>
inline typename codi::TypeTraits<Real>::PassiveReal ceil(const codi::Expression<Real, A>& a) {
  return ceil(a.getValue());
}

}

