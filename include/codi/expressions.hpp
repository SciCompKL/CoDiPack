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

  /*
   * Define the expression interface.
   */

  #include "expressionInterface.hpp"

  /*
   * Now define particular types of expressions, using static
   * polymorphism via the Curiously Recurring Template Pattern
   */

  #include "unaryExpressions.hpp"
  #include "binaryExpressions.hpp"

  /*
   * Implement the operation logic.
   */

  /**
   * @brief Interface for binary elementary operation logic.
   *
   * @tparam Real   The real type used in the active types.
   *
   * Must be implemented for every binary elementary operation.
   * Each implementation is followed by an include of binaryOverloads.tpp with prior #defines.
   *
   * The gradient methods immediately return the jacobie with respect to the first and second argument respectively.
   * The derv methods allow for optimizations during backward traversal of expression trees if active and passive arguments are combined or if the backward paths have computations in common.
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

  /*
   * Now all binary functions are implemented.
   */

  /**
   * @brief Operation logic for f(a,b) = a + b
   */
  template<typename Real>
  struct Add : public BinaryOpInterface<Real> {
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
  #define OPERATION_LOGIC Add
  #define FUNCTION operator +
  #include "binaryOverloads.tpp"

  /**
   * @brief Operation logic for f(a,b) = a - b
   */
  template<typename Real>
  struct Subtract : public BinaryOpInterface<Real> {
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
  #define OPERATION_LOGIC Subtract
  #define FUNCTION operator -
  #include "binaryOverloads.tpp"

  /**
   * @brief Operation logic for f(a,b) = a * b
   */
  template<typename Real>
  struct Multiply : public BinaryOpInterface<Real> {
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
  #define OPERATION_LOGIC Multiply
  #define FUNCTION operator *
  #include "binaryOverloads.tpp"

  /**
   * @brief Operation logic for f(a,b) = a / b
   */
  template<typename Real>
  struct Divide : public BinaryOpInterface<Real> {
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
  #define OPERATION_LOGIC Divide
  #define FUNCTION operator /
  #include "binaryOverloads.tpp"

  using std::atan2;

  /**
   * @brief Operation logic for f(a,b) = atan2(a,b)
   */
  template<typename Real>
  struct Atan2 : public BinaryOpInterface<Real> {
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
  #define OPERATION_LOGIC Atan2
  #define FUNCTION atan2
  #include "binaryOverloads.tpp"

  using std::pow;

  /**
   * @brief Operation logic for f(a,b) = pow(a,b)
   */
  template<typename Real>
  struct Pow : public BinaryOpInterface<Real> {
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
  #define OPERATION_LOGIC Pow
  #define FUNCTION pow
  #include "binaryOverloads.tpp"

  using std::min;

  /**
   * @brief Operation logic for f(a,b) = min(a,b)
   */
  template<typename Real>
  struct Min : public BinaryOpInterface<Real> {
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
  #define OPERATION_LOGIC Min
  #define FUNCTION min
  #include "binaryOverloads.tpp"

  /*
   * Forward of fmin to min
   */
  /**
   * @brief Overload for fmin with the CoDiPack expressions.
   *
   * @param[in] a  The first argument of the operation.
   * @param[in] b  The second argument of the operation.
   *
   * @return BinaryOp11 instanciated for operation logic Min.
   *
   * @tparam Real  The real type used in the active types.
   * @tparam    A  The expression for the first argument of the function.
   * @tparam    B  The expression for the second argument of the function.
   */
  template <typename Real, class A, class B>
  CODI_INLINE BinaryOp11<Real, A, B, Min> fmin(const Expression<Real, A>& a, const Expression<Real, B>& b) {
    return BinaryOp11<Real, A, B, Min>(a.cast(), b.cast());
  }
  /**
   * @brief Overload for fmin with the CoDiPack expressions.
   *
   * @param[in] a  The first argument of the operation.
   * @param[in] b  The second argument of the operation.
   *
   * @return BinaryOp10 instanciated for operation logic Min.
   *
   * @tparam Real  The real type used in the active types.
   * @tparam    A  The expression for the first argument of the function.
   */
  template <typename Real, class A>
  CODI_INLINE BinaryOp10<Real, A, Min> fmin(const Expression<Real, A>& a, const typename TypeTraits<Real>::PassiveReal& b) {
    return BinaryOp10<Real, A, Min>(a.cast(), b);
  }
  /**
   * @brief Overload for fmin with the CoDiPack expressions.
   *
   * @param[in] a  The first argument of the operation.
   * @param[in] b  The second argument of the operation.
   *
   * @return BinaryOp01 instanciated for operation logic Min.
   *
   * @tparam Real  The real type used in the active types.
   * @tparam    B  The expression for the second argument of the function.
   */
  template <typename Real, class B>
  CODI_INLINE BinaryOp01<Real, B, Min> fmin(const typename TypeTraits<Real>::PassiveReal& a, const Expression<Real, B>& b) {
    return BinaryOp01<Real, B, Min>(a, b.cast());
  }

  using std::max;

  /**
   * @brief Operation logic for f(a,b) = max(a,b)
   */
  template<typename Real>
  struct Max : public BinaryOpInterface<Real> {
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
  #define OPERATION_LOGIC Max
  #define FUNCTION max
  #include "binaryOverloads.tpp"

  /*
   * Forward of fmax to max
   */
  /**
   * @brief Overload for fmax with the CoDiPack expressions.
   *
   * @param[in] a  The first argument of the operation.
   * @param[in] b  The second argument of the operation.
   *
   * @return BinaryOp11 instanciated for operation logic Max.
   *
   * @tparam Real  The real type used in the active types.
   * @tparam    A  The expression for the first argument of the function.
   * @tparam    B  The expression for the second argument of the function.
   */
  template <typename Real, class A, class B>
  CODI_INLINE BinaryOp11<Real, A, B, Max> fmax(const Expression<Real, A>& a, const Expression<Real, B>& b) {
    return BinaryOp11<Real, A, B, Max>(a.cast(), b.cast());
  }

  /**
   * @brief Overload for fmax with the CoDiPack expressions.
   *
   * @param[in] a  The first argument of the operation.
   * @param[in] b  The second argument of the operation.
   *
   * @return BinaryOp10 instanciated for operation logic Max.
   *
   * @tparam Real  The real type used in the active types.
   * @tparam    A  The expression for the first argument of the function.
   */
  template <typename Real, class A>
  CODI_INLINE BinaryOp10<Real, A, Max> fmax(const Expression<Real, A>& a, const typename TypeTraits<Real>::PassiveReal& b) {
    return BinaryOp10<Real, A, Max>(a.cast(), b);
  }
  /**
   * @brief Overload for fmax with the CoDiPack expressions.
   *
   * @param[in] a  The first argument of the operation.
   * @param[in] b  The second argument of the operation.
   *
   * @return BinaryOp01 instanciated for operation logic Max.
   *
   * @tparam Real  The real type used in the active types.
   * @tparam    B  The expression for the second argument of the function.
   */
  template <typename Real, class B>
  CODI_INLINE BinaryOp01<Real, B, Max> fmax(const typename TypeTraits<Real>::PassiveReal& a, const Expression<Real, B>& b) {
    return BinaryOp01<Real, B, Max>(a, b.cast());
  }

  using std::copysign;

  /**
   * @brief Operation logic for f(a,b) = copysign(a,b)
   */
  template<typename Real>
  struct Copysign : public BinaryOpInterface<Real> {
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
  #define OPERATION_LOGIC Copysign
  #define FUNCTION copysign
  #include "binaryOverloads.tpp"

  /**
   * @brief Overload for copysignf with the CoDiPack expressions.
   *
   * @param[in] a  The first argument of the operation.
   * @param[in] b  The second argument of the operation.
   *
   * @return BinaryOp11 instanciated for operation logic Copysign.
   *
   * @tparam Real  The real type used in the active types.
   * @tparam    A  The expression for the first argument of the function.
   * @tparam    B  The expression for the second argument of the function.
   */
  template <typename Real, class A, class B>
  CODI_INLINE BinaryOp11<Real, A, B, Copysign> copysignf(const Expression<Real, A>& a, const Expression<Real, B>& b) {
    return BinaryOp11<Real, A, B, Copysign>(a.cast(), b.cast());
  }

  /**
   * @brief Overload for copysignf with the CoDiPack expressions.
   *
   * @param[in] a  The first argument of the operation.
   * @param[in] b  The second argument of the operation.
   *
   * @return BinaryOp10 instanciated for operation logic Copysign.
   *
   * @tparam Real  The real type used in the active types.
   * @tparam    A  The expression for the first argument of the function.
   */
  template <typename Real, class A>
  CODI_INLINE BinaryOp10<Real, A, Copysign> copysignf(const Expression<Real, A>& a, const typename TypeTraits<Real>::PassiveReal& b) {
    return BinaryOp10<Real, A, Copysign>(a.cast(), b);
  }

  /**
   * @brief Overload for copysignf with the CoDiPack expressions.
   *
   * @param[in] a  The first argument of the operation.
   * @param[in] b  The second argument of the operation.
   *
   * @return BinaryOp01 instanciated for operation logic Copysign.
   *
   * @tparam Real  The real type used in the active types.
   * @tparam    B  The expression for the second argument of the function.
   */
  template <typename Real, class B>
  CODI_INLINE BinaryOp01<Real, B, Copysign> copysignf(const typename TypeTraits<Real>::PassiveReal& a, const Expression<Real, B>& b) {
    return BinaryOp01<Real, B, Copysign>(a, b.cast());
  }

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

  /**
   * @brief Operation logic for f(a) = -a
   */
  template<typename Real>
  struct UnaryMinus : public UnaryOpInterface<Real> {
    static CODI_INLINE Real primal(const Real& a) {
      return -a;
    }

    static CODI_INLINE typename TypeTraits<Real>::PassiveReal gradient(const Real& a, const Real& result) {
      CODI_UNUSED(a);
      CODI_UNUSED(result);
      return -1.0;
    }
  };
  #define OPERATION_LOGIC UnaryMinus
  #define FUNCTION operator -
  #include "unaryOverloads.tpp"

  using std::sqrt;

  /**
   * @brief Operation logic for f(a) = sqrt(a)
   */
  template<typename Real>
  struct Sqrt : public UnaryOpInterface<Real> {
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
  #define OPERATION_LOGIC Sqrt
  #define FUNCTION sqrt
  #include "unaryOverloads.tpp"

  using std::cbrt;

  /**
   * @brief Operation logic for f(a) = cbrt(a)
   */
  template<typename Real>
  struct Cbrt : public UnaryOpInterface<Real> {
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
  #define OPERATION_LOGIC Cbrt
  #define FUNCTION cbrt
  #include "unaryOverloads.tpp"

  using std::tanh;

  /**
   * @brief Operation logic for f(a) = tanh(a)
   */
  template<typename Real>
  struct Tanh : public UnaryOpInterface<Real> {
    static CODI_INLINE Real primal(const Real& a) {
      return tanh(a);
    }

    static CODI_INLINE Real gradient(const Real& a, const Real& result) {
      CODI_UNUSED(a);
      return 1 - result * result;
    }
  };
  #define OPERATION_LOGIC Tanh
  #define FUNCTION tanh
  #include "unaryOverloads.tpp"

  using std::log;

  /**
   * @brief Operation logic for f(a) = log(a)
   */
  template<typename Real>
  struct Log : public UnaryOpInterface<Real> {
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
  #define OPERATION_LOGIC Log
  #define FUNCTION log
  #include "unaryOverloads.tpp"

  using std::log10;

  /**
   * @brief Operation logic for f(a) = log10(a)
   */
  template<typename Real>
  struct Log10 : public UnaryOpInterface<Real> {
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
  #define OPERATION_LOGIC Log10
  #define FUNCTION log10
  #include "unaryOverloads.tpp"

  using std::sin;
  using std::cos;

  /**
   * @brief Operation logic for f(a) = sin(a)
   */
  template<typename Real>
  struct Sin : public UnaryOpInterface<Real> {
    static CODI_INLINE Real primal(const Real& a) {
      return sin(a);
    }

    static CODI_INLINE Real gradient(const Real& a, const Real& result) {
      CODI_UNUSED(result);
      return cos(a);
    }
  };
  #define OPERATION_LOGIC Sin
  #define FUNCTION sin
  #include "unaryOverloads.tpp"

  /**
   * @brief Operation logic for f(a) = cos(a)
   */
  template<typename Real>
  struct Cos : public UnaryOpInterface<Real> {
    static CODI_INLINE Real primal(const Real& a) {
      return cos(a);
    }

    static CODI_INLINE Real gradient(const Real& a, const Real& result) {
      CODI_UNUSED(result);
      return -sin(a);
    }
  };
  #define OPERATION_LOGIC Cos
  #define FUNCTION cos
  #include "unaryOverloads.tpp"

  using std::asin;

  /**
   * @brief Operation logic for f(a) = asin(a)
   */
  template<typename Real>
  struct Asin : public UnaryOpInterface<Real> {
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
  #define OPERATION_LOGIC Asin
  #define FUNCTION asin
  #include "unaryOverloads.tpp"

  using std::acos;

  /**
   * @brief Operation logic for f(a) = acos(a)
   */
  template<typename Real>
  struct Acos : public UnaryOpInterface<Real> {
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
  #define OPERATION_LOGIC Acos
  #define FUNCTION acos
  #include "unaryOverloads.tpp"

  using std::atan;

  /**
   * @brief Operation logic for f(a) = atan(a)
   */
  template<typename Real>
  struct Atan : public UnaryOpInterface<Real> {
    static CODI_INLINE Real primal(const Real& a) {
      return atan(a);
    }

    static CODI_INLINE Real gradient(const Real& a, const Real& result) {
      CODI_UNUSED(result);
      return 1.0 / (1 + a * a);
    }
  };
  #define OPERATION_LOGIC Atan
  #define FUNCTION atan
  #include "unaryOverloads.tpp"

  using std::sinh;
  using std::cosh;

  /**
   * @brief Operation logic for f(a) = sinh(a)
   */
  template<typename Real>
  struct Sinh : public UnaryOpInterface<Real> {
    static CODI_INLINE Real primal(const Real& a) {
      return sinh(a);
    }

    static CODI_INLINE Real gradient(const Real& a, const Real& result) {
      CODI_UNUSED(result);
      return cosh(a);
    }
  };
  #define OPERATION_LOGIC Sinh
  #define FUNCTION sinh
  #include "unaryOverloads.tpp"

  /**
   * @brief Operation logic for f(a) = cosh(a)
   */
  template<typename Real>
  struct Cosh : public UnaryOpInterface<Real> {
    static CODI_INLINE Real primal(const Real& a) {
      return cosh(a);
    }

    static CODI_INLINE Real gradient(const Real& a, const Real& result) {
      CODI_UNUSED(result);
      return sinh(a);
    }
  };
  #define OPERATION_LOGIC Cosh
  #define FUNCTION cosh
  #include "unaryOverloads.tpp"

  using std::exp;

  /**
   * @brief Operation logic for f(a) = exp(a)
   */
  template<typename Real>
  struct Exp : public UnaryOpInterface<Real> {
    static CODI_INLINE Real primal(const Real& a) {
      return exp(a);
    }

    static CODI_INLINE Real gradient(const Real& a, const Real& result) {
      CODI_UNUSED(a);
      return result;
    }
  };
  #define OPERATION_LOGIC Exp
  #define FUNCTION exp
  #include "unaryOverloads.tpp"

  using std::atanh;

  /**
   * @brief Operation logic for f(a) = atanh(a)
   */
  template<typename Real>
  struct Atanh : public UnaryOpInterface<Real> {
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
  #define OPERATION_LOGIC Atanh
  #define FUNCTION atanh
  #include "unaryOverloads.tpp"

  using std::abs;

  /**
   * @brief Operation logic for f(a) = abs(a)
   */
  template<typename Real>
  struct Abs : public UnaryOpInterface<Real> {
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
  #define OPERATION_LOGIC Abs
  #define FUNCTION abs
  #include "unaryOverloads.tpp"

  using std::tan;

  /**
   * @brief Operation logic for f(a) = tan(a)
   */
  template<typename Real>
  struct Tan : public UnaryOpInterface<Real> {
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
  #define OPERATION_LOGIC Tan
  #define FUNCTION tan
  #include "unaryOverloads.tpp"

  using std::erf;

  /**
   * @brief Operation logic for f(a) = erf(a)
   */
  template<typename Real>
  struct Erf : public UnaryOpInterface<Real> {
    static CODI_INLINE Real primal(const Real& a) {
      return erf(a);
    }

    static CODI_INLINE Real gradient(const Real& a, const Real& result) {
      CODI_UNUSED(result);
      return 1.128379167095513 * exp( -(a * a) ); // erf'(a) = 2.0 / sqrt(pi) * exp(-a^2)
    }
  };
  #define OPERATION_LOGIC Erf
  #define FUNCTION erf
  #include "unaryOverloads.tpp"

  using std::erfc;

  /**
   * @brief Operation logic for f(a) = erfc(a)
   */
  template<typename Real>
  struct Erfc : public UnaryOpInterface<Real> {
    static CODI_INLINE Real primal(const Real& a) {
      return erfc(a);
    }

    static CODI_INLINE Real gradient(const Real& a, const Real& result) {
      CODI_UNUSED(result);
      return -1.128379167095513 * exp( -(a * a) ); // erfc'(a) = - 2.0 / sqrt(pi) * exp(-a^2)
    }
  };
  #define OPERATION_LOGIC Erfc
  #define FUNCTION erfc
  #include "unaryOverloads.tpp"

  using std::tgamma;

  /**
   * @brief Operation logic for f(a) = tgamma(a)
   */
  template<typename Real>
  struct Tgamma : public UnaryOpInterface<Real> {
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
  #define OPERATION_LOGIC Tgamma
  #define FUNCTION tgamma
  #include "unaryOverloads.tpp"

  /**
   * @brief The fabs function is redirected to abs.
   *
   * @param[in] a The argument of the operation.
   *
   * @return UnaryOp instanciated for operation logic Abs.
   *
   * @tparam Real The real type used in the active types.
   * @tparam A The expression for the argument of the function
   */
  template<typename Real, class A>
  CODI_INLINE UnaryOp<Real, A, Abs> fabs(const codi::Expression<Real, A>& a) {
    return UnaryOp<Real, A, Abs>(a.cast());
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
