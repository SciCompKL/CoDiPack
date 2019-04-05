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
   * @brief Interface for binary elementary operation logic.
   *
   * @tparam Real   The real type used in the active types.
   *
   * Must be implemented for every binary elementary operation.
   * Each implementation is followed by an include of binaryExpression.tpp with prior #defines.
   *
   * The gradient methods immediately return the jacobie with respect to the first and
   * second argument respectively.
   * The derv methods allow for optimizations during backward traversal of expression trees if active and passive arguments are combined.
   *
   * We use the naming scheme dervBB[M].
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
  template<typename Real>
  struct BinaryOpInterface {
    /**
     * @brief primal Primal function call.
     * @param[in] a   First argument of the operation.
     * @param[in] b   Second argument of the operation.
     * @return        The result of the operation.
     */
    static CODI_INLINE Real primal(const Real& a, const Real& b);

    /**
     * @brief gradientA Gradient of operation with respect to first parameter.
     * @tparam ReturnType Indicates that the return type depends on the specific operation.
     * @tparam A      Type of first argument. Kept variable because it can be an expression.
     * @tparam B      Type of second argument. Kept variable because it can be an expression.
     * @param[in] a   First argument of the operation.
     * @param[in] b   Second argument of the operation.
     * @return        Gradient with respect to first parameter.
     */
    template<typename ReturnType, typename A, typename B>
    static CODI_INLINE ReturnType gradientA(const A& a, const B& b, const Real& result);

    /**
     * @brief gradientA Gradient of operation with respect to second parameter.
     * @tparam ReturnType Indicates that the return type depends on the specific operation.
     * @tparam A      Type of first argument. Kept variable because it can be an expression.
     * @tparam B      Type of second argument. Kept variable because it can be an expression.
     * @param[in] a   First argument of the operation.
     * @param[in] b   Second argument of the operation.
     * @return        Gradient with respect to second parameter.
     */
    template<typename ReturnType, typename A, typename B>
    static CODI_INLINE ReturnType gradientB(const A& a, const B& b, const Real& result);

    template<typename Data, typename A, typename B>
    static CODI_INLINE void derv11(Data& data, const A& a, const B& b, const Real& result);

    template<typename Data, typename A, typename B>
    static CODI_INLINE void derv11M(Data& data, const A& a, const B& b, const Real& result, const Real& multiplier);

    template<typename Data, typename A, typename B>
    static CODI_INLINE void derv10(Data& data, const A& a, const typename TypeTraits<Real>::PassiveReal& b, const Real& result);

    template<typename Data, typename A, typename B>
    static CODI_INLINE void derv10M(Data& data, const A& a, const typename TypeTraits<Real>::PassiveReal& b, const Real& result, const Real& multiplier);

    template<typename Data, typename A, typename B>
    static CODI_INLINE void derv01(Data& data, const typename TypeTraits<Real>::PassiveReal& a, const B& b, const Real& result);

    template<typename Data, typename A, typename B>
    static CODI_INLINE void derv01M(Data& data, const typename TypeTraits<Real>::PassiveReal& a, const B& b, const Real& result, const Real& multiplier);
  };

  /**
   * @brief Operator logic for f(a,b) = a + b
   */
  template<typename Real>
  struct AddImpl : public BinaryOpInterface<Real> {
    static CODI_INLINE Real primal(const Real& a, const Real& b) {
      return a + b;
    }

    template<typename A, typename B>
    static CODI_INLINE const typename TypeTraits<Real>::PassiveReal gradientA(const A& a, const B& b, const Real& result) {
      CODI_UNUSED(a);
      CODI_UNUSED(b);
      CODI_UNUSED(result);
      return 1.0;
    }

    template<typename A, typename B>
    static CODI_INLINE  const typename TypeTraits<Real>::PassiveReal gradientB(const A& a, const B& b, const Real& result) {
      CODI_UNUSED(a);
      CODI_UNUSED(b);
      CODI_UNUSED(result);
      return 1.0;
    }

    template<typename Data, typename A, typename B>
    static CODI_INLINE void derv11(Data& data, const A& a, const B& b, const Real& result) {
      CODI_UNUSED(result);
      a.calcGradient(data);
      b.calcGradient(data);
    }

    template<typename Data, typename A, typename B>
    static CODI_INLINE void derv11M(Data& data, const A& a, const B& b, const Real& result, const Real& multiplier) {
      CODI_UNUSED(result);
      a.calcGradient(data, multiplier);
      b.calcGradient(data, multiplier);
    }

    template<typename Data, typename A>
    static CODI_INLINE void derv10(Data& data, const A& a, const typename TypeTraits<Real>::PassiveReal& b, const Real& result) {
      CODI_UNUSED(b);
      CODI_UNUSED(result);
      a.calcGradient(data);
    }

    template<typename Data, typename A>
    static CODI_INLINE void derv10M(Data& data, const A& a, const typename TypeTraits<Real>::PassiveReal& b, const Real& result, const Real& multiplier) {
      CODI_UNUSED(b);
      CODI_UNUSED(result);
      a.calcGradient(data, multiplier);
    }

    template<typename Data, typename B>
    static CODI_INLINE void derv01(Data& data, const typename TypeTraits<Real>::PassiveReal& a, const B& b, const Real& result) {
      CODI_UNUSED(a);
      CODI_UNUSED(result);
      b.calcGradient(data);
    }

    template<typename Data, typename B>
    static CODI_INLINE void derv01M(Data& data, const typename TypeTraits<Real>::PassiveReal& a, const B& b, const Real& result, const Real& multiplier) {
      CODI_UNUSED(a);
      CODI_UNUSED(result);
      b.calcGradient(data, multiplier);
    }
  };
  #define NAME Add
  #define FUNCTION operator +
  #include "binaryExpression.tpp"

   /**
   * @brief Operator logic for f(a,b) = a - b
   */
  template<typename Real>
  struct SubtractImpl : public BinaryOpInterface<Real> {
    static CODI_INLINE Real primal(const Real& a, const Real& b) {
      return a - b;
    }

    template<typename A, typename B>
    static CODI_INLINE  const typename TypeTraits<Real>::PassiveReal gradientA(const A& a, const B& b, const Real& result) {
      CODI_UNUSED(a);
      CODI_UNUSED(b);
      CODI_UNUSED(result);
      return 1.0;
    }

    template<typename A, typename B>
    static CODI_INLINE  const typename TypeTraits<Real>::PassiveReal gradientB(const A& a, const B& b, const Real& result) {
      CODI_UNUSED(a);
      CODI_UNUSED(b);
      CODI_UNUSED(result);
      return -1.0;
    }

    template<typename Data, typename A, typename B>
    static CODI_INLINE void derv11(Data& data, const A& a, const B& b, const Real& result) {
      CODI_UNUSED(result);
      a.calcGradient(data);
      b.calcGradient(data, -1.0);
    }

    template<typename Data, typename A, typename B>
    static CODI_INLINE void derv11M(Data& data, const A& a, const B& b, const Real& result, const Real& multiplier) {
      CODI_UNUSED(result);
      a.calcGradient(data, multiplier);
      b.calcGradient(data, -multiplier);
    }

    template<typename Data, typename A>
    static CODI_INLINE void derv10(Data& data, const A& a, const typename TypeTraits<Real>::PassiveReal& b, const Real& result) {
      CODI_UNUSED(b);
      CODI_UNUSED(result);
      a.calcGradient(data);
    }

    template<typename Data, typename A>
    static CODI_INLINE void derv10M(Data& data, const A& a, const typename TypeTraits<Real>::PassiveReal& b, const Real& result, const Real& multiplier) {
      CODI_UNUSED(b);
      CODI_UNUSED(result);
      a.calcGradient(data, multiplier);
    }

    template<typename Data, typename B>
    static CODI_INLINE void derv01(Data& data, const typename TypeTraits<Real>::PassiveReal& a, const B& b, const Real& result) {
      CODI_UNUSED(a);
      CODI_UNUSED(result);
      b.calcGradient(data, -1.0);
    }

    template<typename Data, typename B>
    static CODI_INLINE void derv01M(Data& data, const typename TypeTraits<Real>::PassiveReal& a, const B& b, const Real& result, const Real& multiplier) {
      CODI_UNUSED(a);
      CODI_UNUSED(result);
      b.calcGradient(data, -multiplier);
    }
  };
  #define NAME Subtract
  #define FUNCTION operator -
  #include "binaryExpression.tpp"

   /**
   * @brief Operator logic for f(a,b) = a * b
   */
  template<typename Real>
  struct MultiplyImpl : public BinaryOpInterface<Real> {
    static CODI_INLINE Real primal(const Real& a, const Real& b) {
      return a * b;
    }

    template<typename A, typename B>
    static CODI_INLINE  const B& gradientA(const A& a, const B& b, const Real& result) {
      CODI_UNUSED(a);
      CODI_UNUSED(result);
      return b;
    }

    template<typename A, typename B>
    static CODI_INLINE  const A& gradientB(const A& a, const B& b, const Real& result) {
      CODI_UNUSED(b);
      CODI_UNUSED(result);
      return a;
    }

    template<typename Data, typename A, typename B>
    static CODI_INLINE void derv11(Data& data, const A& a, const B& b, const Real& result) {
      CODI_UNUSED(result);
      a.calcGradient(data, b.getValue());
      b.calcGradient(data, a.getValue());
    }

    template<typename Data, typename A, typename B>
    static CODI_INLINE void derv11M(Data& data, const A& a, const B& b, const Real& result, const Real& multiplier) {
      CODI_UNUSED(result);
      a.calcGradient(data, b.getValue() * multiplier);
      b.calcGradient(data, a.getValue() * multiplier);
    }

    template<typename Data, typename A>
    static CODI_INLINE void derv10(Data& data, const A& a, const typename TypeTraits<Real>::PassiveReal& b, const Real& result) {
      CODI_UNUSED(result);
      a.calcGradient(data, b);
    }

    template<typename Data, typename A>
    static CODI_INLINE void derv10M(Data& data, const A& a, const typename TypeTraits<Real>::PassiveReal& b, const Real& result, const Real& multiplier) {
      CODI_UNUSED(result);
      a.calcGradient(data, b * multiplier);
    }

    template<typename Data, typename B>
    static CODI_INLINE void derv01(Data& data, const typename TypeTraits<Real>::PassiveReal& a, const B& b, const Real& result) {
      CODI_UNUSED(result);
      b.calcGradient(data, a);
    }

    template<typename Data, typename B>
    static CODI_INLINE void derv01M(Data& data, const typename TypeTraits<Real>::PassiveReal& a, const B& b, const Real& result, const Real& multiplier) {
      CODI_UNUSED(result);
      b.calcGradient(data, a * multiplier);
    }
  };
  #define NAME Multiply
  #define FUNCTION operator *
  #include "binaryExpression.tpp"

   /**
   * @brief Operator logic for f(a,b) = a / b
   */
  template<typename Real>
  struct DivideImpl : public BinaryOpInterface<Real> {
    private:
      /**
       * @brief Helper function which checks if the divisor is zero.
       *
       * Depending of the global option CheckExpressionArguments the function will call
       * CODI_EXCEPTION if b is equal to zero.
       *
       * @tparam Real The real type used in the active type.
       */
      template<typename B>
      static CODI_INLINE void checkArguments(const B& b) {
        if(CheckExpressionArguments) {
          if( 0.0 == TypeTraits<B>::getBaseValue(b)) {
            CODI_EXCEPTION("Division called with divisor of zero.");
          }
        }
      }

    public:
      static CODI_INLINE Real primal(const Real& a, const Real& b) {
        return a / b;
      }

      template<typename A, typename B>
      static CODI_INLINE Real gradientA(const A& a, const B& b, const Real& result) {
        checkArguments(b);
        CODI_UNUSED(a);
        CODI_UNUSED(result);
        return 1.0 / b;
      }

      template<typename A, typename B>
      static CODI_INLINE Real gradientB(const A& a, const B& b, const Real& result) {
        checkArguments(b);
        CODI_UNUSED(a);
        return -result / b;
      }

      template<typename Data, typename A, typename B>
      static CODI_INLINE void derv11(Data& data, const A& a, const B& b, const Real& result) {
        checkArguments(b);
        Real one_over_b = 1.0 / b.getValue();
        a.calcGradient(data, one_over_b);
        b.calcGradient(data, -result * one_over_b);
      }

      template<typename Data, typename A, typename B>
      static CODI_INLINE void derv11M(Data& data, const A& a, const B& b, const Real& result, const Real& multiplier) {
        checkArguments(b);
        Real one_over_b = multiplier / b.getValue();
        a.calcGradient(data, one_over_b);
        b.calcGradient(data, -result * one_over_b);
      }

      template<typename Data, typename A>
      static CODI_INLINE void derv10(Data& data, const A& a, const typename TypeTraits<Real>::PassiveReal& b, const Real& result) {
        CODI_UNUSED(result);
        checkArguments(b);
        Real one_over_b = 1.0 / b;
        a.calcGradient(data, one_over_b);
      }

      template<typename Data, typename A>
      static CODI_INLINE void derv10M(Data& data, const A& a, const typename TypeTraits<Real>::PassiveReal& b, const Real& result, const Real& multiplier) {
        CODI_UNUSED(result);
        checkArguments(b);
        Real one_over_b = multiplier / b;
        a.calcGradient(data, one_over_b);
      }

      template<typename Data, typename B>
      static CODI_INLINE void derv01(Data& data, const typename TypeTraits<Real>::PassiveReal& a, const B& b, const Real& result) {
        CODI_UNUSED(a);
        checkArguments(b);
        Real one_over_b = 1.0 / b.getValue();
        b.calcGradient(data, -result * one_over_b);
      }

      template<typename Data, typename B>
      static CODI_INLINE void derv01M(Data& data, const typename TypeTraits<Real>::PassiveReal& a, const B& b, const Real& result, const Real& multiplier) {
        CODI_UNUSED(a);
        checkArguments(b);
        Real one_over_b = multiplier / b.getValue();
        b.calcGradient(data, -result * one_over_b);
      }
  };
  #define NAME Divide
  #define FUNCTION operator /
  #include "binaryExpression.tpp"

   /**
   * @brief Operator logic for f(a,b) = atan2(a,b)
   */
  using std::atan2;
  template<typename Real>
  struct Atan2Impl : public BinaryOpInterface<Real> {
    private:
      /**
       * @brief Helper function which checks if both arguments are zero.
       *
       * Depending of the global option CheckExpressionArguments the function will call
       * CODI_EXCEPTION if a and b are equal to zero.
       *
       * @tparam A Type of the first argument.
       * @tparam B Type of the second argument.
       */
      template<typename A, typename B>
      static CODI_INLINE void checkArguments(const A& a, const B& b) {
        if(CheckExpressionArguments) {
          if( 0.0 == TypeTraits<A>::getBaseValue(a) &&
              0.0 == TypeTraits<B>::getBaseValue(b)) {
            CODI_EXCEPTION("atan2 called at point (0,0).");
          }
        }
      }

    public:
      static CODI_INLINE Real primal(const Real& a, const Real& b) {
        return atan2(a, b);
      }

      template<typename A, typename B>
      static CODI_INLINE Real gradientA(const A& a, const B& b, const Real& result) {
        CODI_UNUSED(result);
        checkArguments(a, b);
        Real divisor = a * a + b * b;
        divisor = 1.0 / divisor;
        return b * divisor;
      }

      template<typename A, typename B>
      static CODI_INLINE Real gradientB(const A& a, const B& b, const Real& result) {
        CODI_UNUSED(result);
        checkArguments(a, b);
        Real divisor = a * a + b * b;
        divisor = 1.0 / divisor;
        return -a * divisor;
      }

      template<typename Data, typename A, typename B>
      static CODI_INLINE void derv11(Data& data, const A& a, const B& b, const Real& result) {
        CODI_UNUSED(result);
        checkArguments(a, b);
        Real divisor = a.getValue() * a.getValue() + b.getValue() * b.getValue();
        divisor = 1.0 / divisor;
        a.calcGradient(data, b.getValue() * divisor);
        b.calcGradient(data, -a.getValue() * divisor);
      }

      template<typename Data, typename A, typename B>
      static CODI_INLINE void derv11M(Data& data, const A& a, const B& b, const Real& result, const Real& multiplier) {
        CODI_UNUSED(result);
        checkArguments(a, b);
        Real divisor = a.getValue() * a.getValue() + b.getValue() * b.getValue();
        divisor = 1.0 / divisor;
        a.calcGradient(data, multiplier * b.getValue() * divisor);
        b.calcGradient(data, multiplier * -a.getValue() * divisor);
      }

      template<typename Data, typename A>
      static CODI_INLINE void derv10(Data& data, const A& a, const typename TypeTraits<Real>::PassiveReal& b, const Real& result) {
        CODI_UNUSED(result);
        checkArguments(a, b);
        Real divisor = a.getValue() * a.getValue() + b * b;
        divisor = 1.0 / divisor;
        a.calcGradient(data, b * divisor);
      }

      template<typename Data, typename A>
      static CODI_INLINE void derv10M(Data& data, const A& a, const typename TypeTraits<Real>::PassiveReal& b, const Real& result, const Real& multiplier) {
        CODI_UNUSED(result);
        checkArguments(a, b);
        Real divisor = a.getValue() * a.getValue() + b * b;
        divisor = 1.0 / divisor;
        a.calcGradient(data, multiplier * b * divisor);
      }

      template<typename Data, typename B>
      static CODI_INLINE void derv01(Data& data, const typename TypeTraits<Real>::PassiveReal& a, const B& b, const Real& result) {
        CODI_UNUSED(result);
        checkArguments(a, b);
        Real divisor = a * a + b.getValue() * b.getValue();
        divisor = 1.0 / divisor;
        b.calcGradient(data, -a * divisor);
      }

      template<typename Data, typename B>
      static CODI_INLINE void derv01M(Data& data, const typename TypeTraits<Real>::PassiveReal& a, const B& b, const Real& result, const Real& multiplier) {
        CODI_UNUSED(result);
        checkArguments(a, b);
        Real divisor = a * a + b.getValue() * b.getValue();
        divisor = 1.0 / divisor;
        b.calcGradient(data, multiplier * -a * divisor);
      }
  };
  #define NAME Atan2
  #define FUNCTION atan2
  #include "binaryExpression.tpp"

  /**
   * @brief Operator logic for f(a,b) = pow(a,b)
   */
  using std::pow;
  template<typename Real>
  struct PowImpl : public BinaryOpInterface<Real> {
    private:
      /**
       * @brief Helper function which checks if the base of the power is negative.
       *
       * Depending of the global option CheckExpressionArguments the function will call
       * CODI_EXCEPTION if a is smaller than zero.
       *
       * @tparam A Type of the first argument.
       * @tparam B Type of the second argument.
       */
      template<typename A>
      static CODI_INLINE void checkArguments(const A& a) {
        if(CheckExpressionArguments) {
          if( TypeTraits<A>::getBaseValue(a) < 0.0) {
            CODI_EXCEPTION("Negative base for active exponent in pow function. (Value: %0.15e)", TypeTraits<A>::getBaseValue(a));
          }
        }
      }

    public:
      static CODI_INLINE Real primal(const Real& a, const Real& b) {
        return pow(a, b);
      }

      template<typename A, typename B>
      static CODI_INLINE Real gradientA(const A& a, const B& b, const Real& result) {
        CODI_UNUSED(result);
        checkArguments(a);
        return b * pow(a, b - 1.0);
      }

      template<typename A, typename B>
      static CODI_INLINE Real gradientB(const A& a, const B& b, const Real& result) {
        CODI_UNUSED(b);
        checkArguments(a);
        if (a > 0.0) {
          return log(a) * result;
        } else {
          return 0.0;
        }
      }

      template<typename Data, typename A, typename B>
      static CODI_INLINE void derv11(Data& data, const A& a, const B& b, const Real& result) {
        checkArguments(a);
        a.calcGradient(data, b.getValue() * pow(a.getValue(), b.getValue() - 1.0));
        if (a.getValue() > 0.0) {
          b.calcGradient(data, log(a.getValue()) * result);
        }
      }

      template<typename Data, typename A, typename B>
      static CODI_INLINE void derv11M(Data& data, const A& a, const B& b, const Real& result, const Real& multiplier) {
        checkArguments(a);
        a.calcGradient(data, multiplier * b.getValue() * pow(a.getValue(), b.getValue() - 1.0));
        if (a.getValue() > 0.0) {
          b.calcGradient(data, multiplier * log(a.getValue()) * result);
        }
      }

      template<typename Data, typename A>
      static CODI_INLINE void derv10(Data& data, const A& a, const typename TypeTraits<Real>::PassiveReal& b, const Real& result) {
        CODI_UNUSED(result);
        a.calcGradient(data, b * pow(a.getValue(), b - 1.0));
      }

      template<typename Data, typename A>
      static CODI_INLINE void derv10M(Data& data, const A& a, const typename TypeTraits<Real>::PassiveReal& b, const Real& result, const Real& multiplier) {
        CODI_UNUSED(result);
        a.calcGradient(data, multiplier * b * pow(a.getValue(), b - 1.0));
      }

      template<typename Data, typename B>
      static CODI_INLINE void derv01(Data& data, const typename TypeTraits<Real>::PassiveReal& a, const B& b, const Real& result) {
        checkArguments(a);
        if (a > 0.0) {
          b.calcGradient(data, log(a) * result);
        }
      }

      template<typename Data, typename B>
      static CODI_INLINE void derv01M(Data& data, const typename TypeTraits<Real>::PassiveReal& a, const B& b, const Real& result, const Real& multiplier) {
        checkArguments(a);
        if (a > 0.0) {
          b.calcGradient(data, multiplier * log(a) * result);
        }
      }
  };
  #define NAME Pow
  #define FUNCTION pow
  #include "binaryExpression.tpp"

  /**
   * @brief Operator logic for f(a,b) = min(a,b)
   */
  using std::min;
  template<typename Real>
  struct MinImpl : public BinaryOpInterface<Real> {
    static CODI_INLINE Real primal(const Real& a, const Real& b) {
      return min(a, b);
    }

    template<typename A, typename B>
    static CODI_INLINE typename TypeTraits<Real>::PassiveReal gradientA(const A& a, const B& b, const Real& result) {
      CODI_UNUSED(result);
      if(a < b) {
        return 1.0;
      } else {
        return 0.0;
      }
    }

    template<typename A, typename B>
    static CODI_INLINE typename TypeTraits<Real>::PassiveReal gradientB(const A& a, const B& b, const Real& result) {
      CODI_UNUSED(result);
      if(a < b) {
        return 0.0;
      } else {
        return 1.0;
      }
    }

    template<typename Data, typename A, typename B>
    static CODI_INLINE void derv11(Data& data, const A& a, const B& b, const Real& result) {
      CODI_UNUSED(result);
      if(a.getValue() < b.getValue()) {
        a.calcGradient(data);
      } else {
        b.calcGradient(data);
      }
    }

    template<typename Data, typename A, typename B>
    static CODI_INLINE void derv11M(Data& data, const A& a, const B& b, const Real& result, const Real& multiplier) {
      CODI_UNUSED(result);
      if(a.getValue() < b.getValue()) {
        a.calcGradient(data, multiplier);
      } else {
        b.calcGradient(data, multiplier);
      }
    }

    template<typename Data, typename A>
    static CODI_INLINE void derv10(Data& data, const A& a, const typename TypeTraits<Real>::PassiveReal& b, const Real& result) {
      CODI_UNUSED(result);
      if(a.getValue() < b) {
        a.calcGradient(data);
      }
    }

    template<typename Data, typename A>
    static CODI_INLINE void derv10M(Data& data, const A& a, const typename TypeTraits<Real>::PassiveReal& b, const Real& result, const Real& multiplier) {
      CODI_UNUSED(result);
      if(a.getValue() < b) {
        a.calcGradient(data, multiplier);
      }
    }

    template<typename Data, typename B>
    static CODI_INLINE void derv01(Data& data, const typename TypeTraits<Real>::PassiveReal& a, const B& b, const Real& result) {
       CODI_UNUSED(result);
      if(a >= b.getValue()) {
        b.calcGradient(data);
      }
    }

    template<typename Data, typename B>
    static CODI_INLINE void derv01M(Data& data, const typename TypeTraits<Real>::PassiveReal& a, const B& b, const Real& result, const Real& multiplier) {
      CODI_UNUSED(result);
      if(a >= b.getValue()) {
        b.calcGradient(data, multiplier);
      }
    }
  };
  #define NAME Min
  #define FUNCTION min
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

  /**
   * @brief Operator logic for f(a,b) = max(a,b)
   */
  using std::max;
  template<typename Real>
  struct MaxImpl : public BinaryOpInterface<Real> {
    static CODI_INLINE Real primal(const Real& a, const Real& b) {
      return max(a, b);
    }

    template<typename A, typename B>
    static CODI_INLINE typename TypeTraits<Real>::PassiveReal gradientA(const A& a, const B& b, const Real& result) {
      CODI_UNUSED(result);
      if(a > b) {
        return 1.0;
      } else {
        return 0.0;
      }
    }

    template<typename A, typename B>
    static CODI_INLINE typename TypeTraits<Real>::PassiveReal gradientB(const A& a, const B& b, const Real& result) {
      CODI_UNUSED(result);
      if(a > b) {
        return 0.0;
      } else {
        return 1.0;
      }
    }

    template<typename Data, typename A, typename B>
    static CODI_INLINE void derv11(Data& data, const A& a, const B& b, const Real& result) {
      CODI_UNUSED(result);
      if(a.getValue() > b.getValue()) {
        a.calcGradient(data);
      } else {
        b.calcGradient(data);
      }
    }

    template<typename Data, typename A, typename B>
    static CODI_INLINE void derv11M(Data& data, const A& a, const B& b, const Real& result, const Real& multiplier) {
      CODI_UNUSED(result);
      if(a.getValue() > b.getValue()) {
        a.calcGradient(data, multiplier);
      } else {
        b.calcGradient(data, multiplier);
      }
    }

    template<typename Data, typename A>
    static CODI_INLINE void derv10(Data& data, const A& a, const typename TypeTraits<Real>::PassiveReal& b, const Real& result) {
      CODI_UNUSED(result);
      if(a.getValue() > b) {
        a.calcGradient(data);
      }
    }

    template<typename Data, typename A>
    static CODI_INLINE void derv10M(Data& data, const A& a, const typename TypeTraits<Real>::PassiveReal& b, const Real& result, const Real& multiplier) {
      CODI_UNUSED(result);
      if(a.getValue() > b) {
        a.calcGradient(data, multiplier);
      }
    }

    template<typename Data, typename B>
    static CODI_INLINE void derv01(Data& data, const typename TypeTraits<Real>::PassiveReal& a, const B& b, const Real& result) {
      CODI_UNUSED(result);
      if(a <= b.getValue()) {
        b.calcGradient(data);
      }
    }

    template<typename Data, typename B>
    static CODI_INLINE void derv01M(Data& data, const typename TypeTraits<Real>::PassiveReal& a, const B& b, const Real& result, const Real& multiplier) {
      CODI_UNUSED(result);
      if(a <= b.getValue()) {
        b.calcGradient(data, multiplier);
      }
    }
  };
  #define NAME Max
  #define FUNCTION max
  #include "binaryExpression.tpp"

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

  /**
   * @brief Operator logic for f(a,b) = copysign(a,b)
   */
  using std::copysign;
  template<typename Real>
  struct CopysignImpl : public BinaryOpInterface<Real> {
    static CODI_INLINE Real primal(const Real& a, const Real& b) {
      return copysign(a, b);
    }

    template<typename A, typename B>
    static CODI_INLINE typename TypeTraits<Real>::PassiveReal gradientA(const A& a, const B& b, const Real& result) {
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

    template<typename A, typename B>
    static CODI_INLINE typename TypeTraits<Real>::PassiveReal gradientB(const A& a, const B& b, const Real& result) {
      CODI_UNUSED(a);
      CODI_UNUSED(b);
      CODI_UNUSED(result);
      return 0.0;
    }

    template<typename Data, typename A, typename B>
    static CODI_INLINE void derv11(Data& data, const A& a, const B& b, const Real& result) {
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

    template<typename Data, typename A, typename B>
    static CODI_INLINE void derv11M(Data& data, const A& a, const B& b, const Real& result, const Real& multiplier) {
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

    template<typename Data, typename A>
    static CODI_INLINE void derv10(Data& data, const A& a, const typename TypeTraits<Real>::PassiveReal& b, const Real& result) {
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

    template<typename Data, typename A>
    static CODI_INLINE void derv10M(Data& data, const A& a, const typename TypeTraits<Real>::PassiveReal& b, const Real& result, const Real& multiplier) {
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

    template<typename Data, typename B>
    static CODI_INLINE void derv01(Data& data, const typename TypeTraits<Real>::PassiveReal& a, const B& b, const Real& result) {
      CODI_UNUSED(a);
      CODI_UNUSED(b);
      CODI_UNUSED(result);
      b.calcGradient(data, 0.0);
    }

    template<typename Data, typename B>
    static CODI_INLINE void derv01M(Data& data, const typename TypeTraits<Real>::PassiveReal& a, const B& b, const Real& result, const Real& multiplier) {
      CODI_UNUSED(a);
      CODI_UNUSED(b);
      CODI_UNUSED(result);
      b.calcGradient(data, 0.0);
    }
  };
  #define NAME Copysign
  #define FUNCTION copysign
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

  /**
   * @brief Interface for unary elementary operation logic.
   *
   * @tparam Real   The real type used in the active types.
   *
   * Must be implemented for every unary elementary operation.
   * Each implementation is followed by an include of unaryExpression.tpp with prior #defines.
   */
  template<typename Real>
  struct UnaryOpInterface {
    /**
     * @brief primal Primal function call.
     * @param[in] a   The argument of the operation.
     * @return        The result of the operation.
     */
    static CODI_INLINE Real primal(const Real& a);

    /**
     * @brief gradient Gradient of the operation.
     * @tparam ReturnType Indicates that the return type depends on the specific elementary operation.
     * @param[in] a       Argument of the operation.
     * @param[in] result  Result of the primal function call.
     * @return            Gradient value.
     */
    template<typename ReturnType>
    static CODI_INLINE ReturnType gradient(const Real& a, const Real& result);
  };

  /*
   * Now all unary functions are implemented.
   */

  template<typename Real>
  struct UnaryMinusImpl : public UnaryOpInterface<Real> {
    static CODI_INLINE Real primal(const Real& a) {
      return -a;
    }

    static CODI_INLINE typename TypeTraits<Real>::PassiveReal gradient(const Real& a, const Real& result) {
      CODI_UNUSED(a);
      CODI_UNUSED(result);
      return -1.0;
    }
  };
  #define NAME UnaryMinus
  #define FUNCTION operator -
  #include "unaryExpression.tpp"

  using std::sqrt;
  template<typename Real>
  struct SqrtImpl : public UnaryOpInterface<Real> {
    static CODI_INLINE Real primal(const Real& a) {
      return sqrt(a);
    }

    static CODI_INLINE Real gradient(const Real& a, const Real& result) {
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
  };
  #define NAME Sqrt
  #define FUNCTION sqrt
  #include "unaryExpression.tpp"

  using std::cbrt;
  template<typename Real>
  struct CbrtImpl : public UnaryOpInterface<Real> {
    static CODI_INLINE Real primal(const Real& a) {
      return cbrt(a);
    }

    static CODI_INLINE Real gradient(const Real& a, const Real& result) {
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
  };
  #define NAME Cbrt
  #define FUNCTION cbrt
  #include "unaryExpression.tpp"

  using std::tanh;
  template<typename Real>
  struct TanhImpl : public UnaryOpInterface<Real> {
    static CODI_INLINE Real primal(const Real& a) {
      return tanh(a);
    }

    static CODI_INLINE Real gradient(const Real& a, const Real& result) {
      CODI_UNUSED(a);
      return 1 - result * result;
    }
  };
  #define NAME Tanh
  #define FUNCTION tanh
  #include "unaryExpression.tpp"

  using std::log;
  template<typename Real>
  struct LogImpl : public UnaryOpInterface<Real> {
    static CODI_INLINE Real primal(const Real& a) {
      return log(a);
    }

    static CODI_INLINE Real gradient(const Real& a, const Real& result) {
      CODI_UNUSED(result);
      if(CheckExpressionArguments) {
        if(0.0 > TypeTraits<Real>::getBaseValue(a)) {
          CODI_EXCEPTION("Logarithm of negative value or zero.(Value: %0.15e)", TypeTraits<Real>::getBaseValue(a));
        }
      }
      return 1.0 / a;
    }
  };
  #define NAME Log
  #define FUNCTION log
  #include "unaryExpression.tpp"

  using std::log10;
  template<typename Real>
  struct Log10Impl : public UnaryOpInterface<Real> {
    static CODI_INLINE Real primal(const Real& a) {
      return log10(a);
    }

    static CODI_INLINE Real gradient(const Real& a, const Real& result) {
      CODI_UNUSED(result);
      if(CheckExpressionArguments) {
        if(0.0 > TypeTraits<Real>::getBaseValue(a)) {
          CODI_EXCEPTION("Logarithm of negative value or zero.(Value: %0.15e)", TypeTraits<Real>::getBaseValue(a));
        }
      }
      return 0.434294481903252 / a;
    }
  };
  #define NAME Log10
  #define FUNCTION log10
  #include "unaryExpression.tpp"

  using std::sin;
  using std::cos;
  template<typename Real>
  struct SinImpl : public UnaryOpInterface<Real> {
    static CODI_INLINE Real primal(const Real& a) {
      return sin(a);
    }

    static CODI_INLINE Real gradient(const Real& a, const Real& result) {
      CODI_UNUSED(result);
      return cos(a);
    }
  };
  #define NAME Sin
  #define FUNCTION sin
  #include "unaryExpression.tpp"

  template<typename Real>
  struct CosImpl : public UnaryOpInterface<Real> {
    static CODI_INLINE Real primal(const Real& a) {
      return cos(a);
    }

    static CODI_INLINE Real gradient(const Real& a, const Real& result) {
      CODI_UNUSED(result);
      return -sin(a);
    }
  };
  #define NAME Cos
  #define FUNCTION cos
  #include "unaryExpression.tpp"

  using std::asin;
  template<typename Real>
  struct AsinImpl : public UnaryOpInterface<Real> {
    static CODI_INLINE Real primal(const Real& a) {
      return asin(a);
    }

    static CODI_INLINE Real gradient(const Real& a, const Real& result) {
      CODI_UNUSED(result);
      if(CheckExpressionArguments) {
        if(TypeTraits<Real>::getBaseValue(a) <= -1.0 || 1.0 <= TypeTraits<Real>::getBaseValue(a)) {
          CODI_EXCEPTION("asin outside of (-1, 1).(Value: %0.15e)", TypeTraits<Real>::getBaseValue(a));
        }
      }
      return 1.0 / sqrt(1.0 - a * a);
    }
  };
  #define NAME Asin
  #define FUNCTION asin
  #include "unaryExpression.tpp"

  using std::acos;
  template<typename Real>
  struct AcosImpl : public UnaryOpInterface<Real> {
    static CODI_INLINE Real primal(const Real& a) {
      return acos(a);
    }

    static CODI_INLINE Real gradient(const Real& a, const Real& result) {
      CODI_UNUSED(result);
      if(CheckExpressionArguments) {
        if(TypeTraits<Real>::getBaseValue(a) <= -1.0 || 1.0 <= TypeTraits<Real>::getBaseValue(a)) {
          CODI_EXCEPTION("acos outside of (-1, 1).(Value: %0.15e)", TypeTraits<Real>::getBaseValue(a));
        }
      }
      return -1.0 / sqrt(1.0 - a * a);
    }
  };
  #define NAME Acos
  #define FUNCTION acos
  #include "unaryExpression.tpp"

  using std::atan;
  template<typename Real>
  struct AtanImpl : public UnaryOpInterface<Real> {
    static CODI_INLINE Real primal(const Real& a) {
      return atan(a);
    }

    static CODI_INLINE Real gradient(const Real& a, const Real& result) {
      CODI_UNUSED(result);
      return 1.0 / (1 + a * a);
    }
  };
  #define NAME Atan
  #define FUNCTION atan
  #include "unaryExpression.tpp"

  using std::sinh;
  using std::cosh;
  template<typename Real>
  struct SinhImpl : public UnaryOpInterface<Real> {
    static CODI_INLINE Real primal(const Real& a) {
      return sinh(a);
    }

    static CODI_INLINE Real gradient(const Real& a, const Real& result) {
      CODI_UNUSED(result);
      return cosh(a);
    }
  };
  #define NAME Sinh
  #define FUNCTION sinh
  #include "unaryExpression.tpp"

  template<typename Real>
  struct CoshImpl : public UnaryOpInterface<Real> {
    static CODI_INLINE Real primal(const Real& a) {
      return cosh(a);
    }

    static CODI_INLINE Real gradient(const Real& a, const Real& result) {
      CODI_UNUSED(result);
      return sinh(a);
    }
  };
  #define NAME Cosh
  #define FUNCTION cosh
  #include "unaryExpression.tpp"

  using std::exp;
  template<typename Real>
  struct ExpImpl : public UnaryOpInterface<Real> {
    static CODI_INLINE Real primal(const Real& a) {
      return exp(a);
    }

    static CODI_INLINE Real gradient(const Real& a, const Real& result) {
      CODI_UNUSED(a);
      return result;
    }
  };
  #define NAME Exp
  #define FUNCTION exp
  #include "unaryExpression.tpp"

  using std::atanh;
  template<typename Real>
  struct AtanhImpl : public UnaryOpInterface<Real> {
    static CODI_INLINE Real primal(const Real& a) {
      return atanh(a);
    }

    static CODI_INLINE Real gradient(const Real& a, const Real& result) {
      CODI_UNUSED(result);
      if(CheckExpressionArguments) {
        if(TypeTraits<Real>::getBaseValue(a) <= -1.0 || 1.0 <= TypeTraits<Real>::getBaseValue(a)) {
          CODI_EXCEPTION("atanh outside of (-1, 1).(Value: %0.15e)", TypeTraits<Real>::getBaseValue(a));
        }
      }
      return 1.0 / (1 - a * a);
    }
  };
  #define NAME Atanh
  #define FUNCTION atanh
  #include "unaryExpression.tpp"

  using std::abs;
  template<typename Real>
  struct AbsImpl : public UnaryOpInterface<Real> {
    static CODI_INLINE Real primal(const Real& a) {
      return abs(a);
    }

    static CODI_INLINE Real gradient(const Real& a, const Real& result) {
      CODI_UNUSED(result);
      if(a < 0.0) {
        return (Real)-1.0;
      } else if(a > 0.0) {
        return (Real)1.0;
      } else {
        return (Real)0.0;
      }
    }
  };
  #define NAME Abs
  #define FUNCTION abs
  #include "unaryExpression.tpp"

  using std::tan;
  template<typename Real>
  struct TanImpl : public UnaryOpInterface<Real> {
    static CODI_INLINE Real primal(const Real& a) {
      return tan(a);
    }

    static CODI_INLINE Real gradient(const Real& a, const Real& result) {
      CODI_UNUSED(result);
      if(CheckExpressionArguments) {
        if(0.0 == cos(TypeTraits<Real>::getBaseValue(a))) {
          CODI_EXCEPTION("Tan evaluated at (0.5  + i) * PI.(Value: %0.15e)", TypeTraits<Real>::getBaseValue(a));
        }
      }
      Real tmp = 1.0 / cos(a);
      return tmp * tmp;
    }
  };
  #define NAME Tan
  #define FUNCTION tan
  #include "unaryExpression.tpp"

  using std::erf;
  template<typename Real>
  struct ErfImpl : public UnaryOpInterface<Real> {
    static CODI_INLINE Real primal(const Real& a) {
      return erf(a);
    }

    static CODI_INLINE Real gradient(const Real& a, const Real& result) {
      CODI_UNUSED(result);
      return 1.128379167095513 * exp( -(a * a) ); // erf'(a) = 2.0 / sqrt(pi) * exp(-a^2)
    }
  };
  #define NAME Erf
  #define FUNCTION erf
  #include "unaryExpression.tpp"

  using std::erfc;
  template<typename Real>
  struct ErfcImpl : public UnaryOpInterface<Real> {
    static CODI_INLINE Real primal(const Real& a) {
      return erfc(a);
    }

    static CODI_INLINE Real gradient(const Real& a, const Real& result) {
      CODI_UNUSED(result);
      return -1.128379167095513 * exp( -(a * a) ); // erfc'(a) = - 2.0 / sqrt(pi) * exp(-a^2)
    }
  };
  #define NAME Erfc
  #define FUNCTION erfc
  #include "unaryExpression.tpp"

  using std::tgamma;
  template<typename Real>
  struct TgammaImpl : public UnaryOpInterface<Real> {
    static CODI_INLINE Real primal(const Real& a) {
      return tgamma(a);
    }

    static CODI_INLINE Real gradient(const Real& a, const Real& result) {
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
  };
  #define NAME Tgamma
  #define FUNCTION tgamma
  #include "unaryExpression.tpp"

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
