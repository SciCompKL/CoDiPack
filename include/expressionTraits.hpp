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

#include <cstddef>

#include "expressions.h"

namespace codi {
  template<class A>
  struct ExpressionTraits {
    /**
     * @brief The maximum number of active variables for the expression
     *
     * This value ca be used to determine how many variables can be active in an expression.
     * For evary expression a specialization has to be defined which assigns a value to
     * the varaible.
     */
    static const size_t maxActiveVariable;
  };

  // Macro for expresions with two arguments
# define ADEPT_DEFINE_BINARY_TRAIT(OP)                        \
    template<typename Real, typename A, typename B>           \
    struct ExpressionTraits<OP<Real, A, B> > {                \
      static const size_t maxActiveVariables =                \
           ExpressionTraits<A>::maxActiveVariables            \
         + ExpressionTraits<B>::maxActiveVariables;           \
    };

  ADEPT_DEFINE_BINARY_TRAIT(Add)
  ADEPT_DEFINE_BINARY_TRAIT(Subtract)
  ADEPT_DEFINE_BINARY_TRAIT(Multiply)
  ADEPT_DEFINE_BINARY_TRAIT(Divide)
  ADEPT_DEFINE_BINARY_TRAIT(Pow)
  ADEPT_DEFINE_BINARY_TRAIT(Atan2)

# undef ADEPT_DEFINE_BINARY_TRAIT

  // Macro for expressions with one argument
# define ADEPT_DEFINE_UNARY_TRAIT(OP)                         \
    template<typename Real, typename A>                       \
    struct ExpressionTraits<OP<Real, A> > {                   \
      static const size_t maxActiveVariables                  \
         = ExpressionTraits<A>::maxActiveVariables;           \
    };

  ADEPT_DEFINE_UNARY_TRAIT(ScalarAdd)
  ADEPT_DEFINE_UNARY_TRAIT(ScalarSubtract)
  ADEPT_DEFINE_UNARY_TRAIT(ScalarMultiply)
  ADEPT_DEFINE_UNARY_TRAIT(ScalarDivide)
  ADEPT_DEFINE_UNARY_TRAIT(UnaryMinus)
  ADEPT_DEFINE_UNARY_TRAIT(Exp)
  ADEPT_DEFINE_UNARY_TRAIT(Tan)
  ADEPT_DEFINE_UNARY_TRAIT(Log)
  ADEPT_DEFINE_UNARY_TRAIT(Log10)
  ADEPT_DEFINE_UNARY_TRAIT(Sqrt)
  ADEPT_DEFINE_UNARY_TRAIT(Sin)
  ADEPT_DEFINE_UNARY_TRAIT(Cos)
  ADEPT_DEFINE_UNARY_TRAIT(Asin)
  ADEPT_DEFINE_UNARY_TRAIT(Acos)
  ADEPT_DEFINE_UNARY_TRAIT(Atan)
  ADEPT_DEFINE_UNARY_TRAIT(Sinh)
  ADEPT_DEFINE_UNARY_TRAIT(Cosh)
  ADEPT_DEFINE_UNARY_TRAIT(Tanh)
  ADEPT_DEFINE_UNARY_TRAIT(Abs)
  ADEPT_DEFINE_UNARY_TRAIT(PowScalarExponent)
  ADEPT_DEFINE_UNARY_TRAIT(PowScalarBase)
  ADEPT_DEFINE_UNARY_TRAIT(Atanh)
  ADEPT_DEFINE_UNARY_TRAIT(Atan2Scalar1)
  ADEPT_DEFINE_UNARY_TRAIT(Atan2Scalar2)

# undef ADEPT_DEFINE_UNARY_TRAIT
}
