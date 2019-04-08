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

  /**
   * @brief Specialization for BinaryOp11. @tparam Real The real type used in the active types.
   *
   * @tparam A The expression for the first argumen of the function.
   * @tparam B The expression for the second argument of the function.
   * @tparam Impl Operation logic implementation.
   */
  template<typename Real, typename A, typename B, template<typename> class Impl>
  struct ExpressionTraits<BinaryOp11<Real, A, B, Impl> > {
    /** @brief Number of maximum active variables is the sum of the active variables from both arguments. */
    static const size_t maxActiveVariables =
         ExpressionTraits<A>::maxActiveVariables
       + ExpressionTraits<B>::maxActiveVariables;
    /** @brief Number of maximum passive variables is the sum of the passive variables from both arguments. */
    static const size_t maxConstantVariables =
         ExpressionTraits<A>::maxConstantVariables
       + ExpressionTraits<B>::maxConstantVariables;
  };

  /**
   * @brief Specialization for BinaryOp10  with only the first argument active.
   *
   * @tparam Real The real type used in the active types.
   * @tparam A The expression for the first argument of the function.
   * @tparam Impl Operation logic implementation.
   */
  template<typename Real, typename A, template<typename> class Impl>
  struct ExpressionTraits<BinaryOp10<Real, A, Impl> > {
    /** @brief Number of maximum active variables is the number of active variables from the first argument. */
    static const size_t maxActiveVariables =
         ExpressionTraits<A>::maxActiveVariables;
    /** @brief Number of maximum passive variables is the number of passive variables from the first argument plus the passive value from this expression. */
    static const size_t maxConstantVariables =
         1 + ExpressionTraits<A>::maxConstantVariables;
  };

  /**
   * @brief Specialization for BinaryOp01 with only the second argument active.
   *
   * @tparam Real The real type used in the active types.
   * @tparam B The expression for the second argument of the function.
   * @tparam Impl Operation logic implementation.
   */
  template<typename Real, typename B, template<typename> class Impl>
  struct ExpressionTraits<BinaryOp01<Real, B, Impl> > {
    /** @brief Number of maximum active variables is the number of active variables from the second argument. */
    static const size_t maxActiveVariables =
         ExpressionTraits<B>::maxActiveVariables;
    /** @brief Number of maximum passive variables is the number of passive variables from the second argument plus the passive value from this expression. */
    static const size_t maxConstantVariables =
         1 + ExpressionTraits<B>::maxConstantVariables;
  };

  /**
   * @brief Specialization for UnaryOp.
   *
   * @tparam Real The real type used in the active types.
   * @tparam A The expression for the argument of the function.
   * @tparam Impl Operation logic implementation.
   */
  template<typename Real, typename A, template<typename> class Impl>
  struct ExpressionTraits<UnaryOp<Real, A, Impl> > {
    /** @brief Number of maximum active variables is the number of active variables from the expression in the argument. */
    static const size_t maxActiveVariables
       = ExpressionTraits<A>::maxActiveVariables;
    /** @brief Number of maximum passive variables is the number of passive variables from the expression in the argument. */
    static const size_t maxConstantVariables
       = ExpressionTraits<A>::maxConstantVariables;
  };
}
