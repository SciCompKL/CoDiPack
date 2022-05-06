/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2021 Chair for Scientific Computing (SciComp), TU Kaiserslautern
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
 * Authors:
 *  - SciComp, TU Kaiserslautern:
 *     Max Sagebaum
 *     Tim Albring
 *     Johannes Blühdorn
 */

#pragma once

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
  };
}
