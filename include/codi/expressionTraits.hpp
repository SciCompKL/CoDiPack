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

#include <cstddef>

#include "expressions.hpp"

/**
 * @brief Global namespace for CoDiPack - Code Differentiation Package
 */
namespace codi {

  /**
   * @brief Information about the expression.
   *
   * The class contains information about specific expressions. This information
   * is used by CoDiPack to run the differentiation.
   *
   * @tparam A The expression for which the information is needed.
   */
  template<class A>
  struct ExpressionTraits {
    /**
     * @brief The maximum number of active variables for the expression
     *
     * This value can be used to determine how many variables can be active in an expression.
     * For every expression a specialization has to be defined which assigns a value to
     * the variable.
     */
    static const size_t maxActiveVariable;
    /**
     * @brief The maximum number of passive variables for the expression
     *
     * This value can be used to determine how many variables are passive in an expression.
     * For every expression a specialization has to be defined which assigns a value to
     * the variable.
     */
    static const size_t maxPassiveVariable;
  };

  // Macro for expressions with two arguments
# define CODI_DEFINE_BINARY_TRAIT(OP)                         \
    /** @brief Specialization for OP. @tparam Real The real type used in the active types. @tparam A The expression for the first argument of the function @tparam B The expression for the second argument of the function*/ \
    template<typename Real, typename A, typename B>           \
    struct ExpressionTraits<OP ## 11<Real, A, B> > {          \
      /** @brief Number of maximum active variables is the sum of the active variables from both arguments. */ \
      static const size_t maxActiveVariables =                \
           ExpressionTraits<A>::maxActiveVariables            \
         + ExpressionTraits<B>::maxActiveVariables;           \
      /** @brief Number of maximum passive variables is the sum of the passive variables from both arguments. */ \
      static const size_t maxConstantVariables =                \
           ExpressionTraits<A>::maxConstantVariables            \
         + ExpressionTraits<B>::maxConstantVariables;           \
    };                                                         \
    /** @brief Specialization for OP  with only the first argument active. @tparam Real The real type used in the active types. @tparam A The expression for the first argument of the function */ \
    template<typename Real, typename A>                       \
    struct ExpressionTraits<OP ## 10<Real, A> > {             \
      /** @brief Number of maximum active variables is the number of active variables from the first argument. */ \
      static const size_t maxActiveVariables =                \
           ExpressionTraits<A>::maxActiveVariables;           \
      /** @brief Number of maximum passive variables is the number of passive variables from the first argument plus the passive value from this expression. */ \
      static const size_t maxConstantVariables =                \
           1 + ExpressionTraits<A>::maxConstantVariables;       \
    };                                                         \
    /** @brief Specialization for OP with only the second argument active. @tparam Real The real type used in the active types. @tparam B The expression for the second argument of the function */ \
    template<typename Real, typename B>                       \
    struct ExpressionTraits<OP ## 01<Real, B> > {             \
      /** @brief Number of maximum active variables is the number of active variables from the second argument. */ \
      static const size_t maxActiveVariables =                \
           ExpressionTraits<B>::maxActiveVariables;           \
      /** @brief Number of maximum passive variables is the number of passive variables from the second argument plus the passive value from this expression. */ \
      static const size_t maxConstantVariables =               \
           1 + ExpressionTraits<B>::maxConstantVariables;      \
    };

  CODI_DEFINE_BINARY_TRAIT(Add)
  CODI_DEFINE_BINARY_TRAIT(Subtract)
  CODI_DEFINE_BINARY_TRAIT(Multiply)
  CODI_DEFINE_BINARY_TRAIT(Divide)
  CODI_DEFINE_BINARY_TRAIT(Pow)
  CODI_DEFINE_BINARY_TRAIT(Atan2)
  CODI_DEFINE_BINARY_TRAIT(Min)
  CODI_DEFINE_BINARY_TRAIT(Max)
  CODI_DEFINE_BINARY_TRAIT(Copysign)

# undef CODI_DEFINE_BINARY_TRAIT

  // Macro for expressions with one argument
# define CODI_DEFINE_UNARY_TRAIT(OP)                          \
    /** @brief Specialization for OP. @tparam Real The real type used in the active types. @tparam A The expression for the argument of the function */ \
    template<typename Real, typename A>                       \
    struct ExpressionTraits<OP<Real, A> > {                   \
      /** @brief Number of maximum active variables is the number of active variables from the expression in the argument. */ \
      static const size_t maxActiveVariables                  \
         = ExpressionTraits<A>::maxActiveVariables;           \
      /** @brief Number of maximum passive variables is the number of passive variables from the expression in the argument. */ \
      static const size_t maxConstantVariables                  \
         = ExpressionTraits<A>::maxConstantVariables;           \
    };

  CODI_DEFINE_UNARY_TRAIT(UnaryMinus)
  CODI_DEFINE_UNARY_TRAIT(Exp)
  CODI_DEFINE_UNARY_TRAIT(Tan)
  CODI_DEFINE_UNARY_TRAIT(Log)
  CODI_DEFINE_UNARY_TRAIT(Log10)
  CODI_DEFINE_UNARY_TRAIT(Sqrt)
  CODI_DEFINE_UNARY_TRAIT(Cbrt)
  CODI_DEFINE_UNARY_TRAIT(Sin)
  CODI_DEFINE_UNARY_TRAIT(Cos)
  CODI_DEFINE_UNARY_TRAIT(Asin)
  CODI_DEFINE_UNARY_TRAIT(Acos)
  CODI_DEFINE_UNARY_TRAIT(Atan)
  CODI_DEFINE_UNARY_TRAIT(Sinh)
  CODI_DEFINE_UNARY_TRAIT(Cosh)
  CODI_DEFINE_UNARY_TRAIT(Tanh)
  CODI_DEFINE_UNARY_TRAIT(Abs)
  CODI_DEFINE_UNARY_TRAIT(Tgamma)
  CODI_DEFINE_UNARY_TRAIT(Atanh)
  CODI_DEFINE_UNARY_TRAIT(Erf)
  CODI_DEFINE_UNARY_TRAIT(Erfc)
# undef CODI_DEFINE_UNARY_TRAIT
}
