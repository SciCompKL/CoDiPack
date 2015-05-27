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
# define CODI_DEFINE_BINARY_TRAIT(OP)                         \
    template<typename Real, typename A, typename B>           \
    struct ExpressionTraits<OP ## 11<Real, A, B> > {          \
      static const size_t maxActiveVariables =                \
           ExpressionTraits<A>::maxActiveVariables            \
         + ExpressionTraits<B>::maxActiveVariables;           \
    };                                                        \
    template<typename Real, typename A>                       \
    struct ExpressionTraits<OP ## 10<Real, A> > {             \
      static const size_t maxActiveVariables =                \
           ExpressionTraits<A>::maxActiveVariables;           \
    };                                                        \
    template<typename Real, typename B>                       \
    struct ExpressionTraits<OP ## 01<Real, B> > {             \
      static const size_t maxActiveVariables =                \
           ExpressionTraits<B>::maxActiveVariables;           \
    };

  CODI_DEFINE_BINARY_TRAIT(Add)
  CODI_DEFINE_BINARY_TRAIT(Subtract)
  CODI_DEFINE_BINARY_TRAIT(Multiply)
  CODI_DEFINE_BINARY_TRAIT(Divide)
  CODI_DEFINE_BINARY_TRAIT(Pow)
  CODI_DEFINE_BINARY_TRAIT(Atan2)
  CODI_DEFINE_BINARY_TRAIT(Min)
  CODI_DEFINE_BINARY_TRAIT(Max)

# undef CODI_DEFINE_BINARY_TRAIT

  // Macro for expressions with one argument
# define CODI_DEFINE_UNARY_TRAIT(OP)                          \
    template<typename Real, typename A>                       \
    struct ExpressionTraits<OP<Real, A> > {                   \
      static const size_t maxActiveVariables                  \
         = ExpressionTraits<A>::maxActiveVariables;           \
    };

  CODI_DEFINE_UNARY_TRAIT(UnaryMinus)
  CODI_DEFINE_UNARY_TRAIT(Exp)
  CODI_DEFINE_UNARY_TRAIT(Tan)
  CODI_DEFINE_UNARY_TRAIT(Log)
  CODI_DEFINE_UNARY_TRAIT(Log10)
  CODI_DEFINE_UNARY_TRAIT(Sqrt)
  CODI_DEFINE_UNARY_TRAIT(Sin)
  CODI_DEFINE_UNARY_TRAIT(Cos)
  CODI_DEFINE_UNARY_TRAIT(Asin)
  CODI_DEFINE_UNARY_TRAIT(Acos)
  CODI_DEFINE_UNARY_TRAIT(Atan)
  CODI_DEFINE_UNARY_TRAIT(Sinh)
  CODI_DEFINE_UNARY_TRAIT(Cosh)
  CODI_DEFINE_UNARY_TRAIT(Tanh)
  CODI_DEFINE_UNARY_TRAIT(Abs)
  CODI_DEFINE_UNARY_TRAIT(Atanh)
# undef CODI_DEFINE_UNARY_TRAIT
}
