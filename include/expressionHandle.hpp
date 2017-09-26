/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2017 Chair for Scientific Computing (SciComp), TU Kaiserslautern
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

#include "typeTraits.hpp"

/**
 * @brief Global namespace for CoDiPack - Code Differentiation Package
 */
namespace codi {

  /**
   * @brief Handle for an expression object.
   *
   * The handle stores information about the expression and the function for the adjoint evaluation.
   *
   * @tparam GradientData  The type of the gradient values that are used in the tape.
   * @tparam         Real  The floating point type of the values that are used in the tape.
   * @tparam    IndexType  The types for the management of the data.
   */
  template<typename GradientData, typename Real, typename IndexType>
  class ExpressionHandle {
    private:
      /**
       * @brief The passive type that corresponds to the real type.
       */
      typedef typename TypeTraits<Real>::PassiveReal PassiveReal;
    public:

      /**
       * @brief The function definition for the reverse evaluation of a statement
       */
      typedef void (*StatementFuncPointer)(const Real& seed, const IndexType* indices, const PassiveReal* passiveValues, const Real* primalValues, Real* adjointValues);

      /**
       * @brief The function pointer to the reverse evaluation function
       */
      const StatementFuncPointer adjointFunc;

      /**
       * @brief The maximum number of active variables in the statement.
       *
       * The number is equal to all the active reals in the statement.
       */
      const size_t maxActiveVariables;

      /**
       * @brief The number of constant values in the statement.
       *
       * The number is equal to all the passive reals in the statement.
       */
      const size_t maxConstantVariables;

      /**
       * @brief Creates the function handle object
       *
       * @param[in]          adjointFunc  The function pointer for the adjoint evaluation.
       * @param[in]   maxActiveVariables  The number of active variables in the statement.
       * @param[in] maxConstantVariables  The number of constant variables in the statement.
       */
      ExpressionHandle(const StatementFuncPointer adjointFunc, const size_t maxActiveVariables, const size_t maxConstantVariables) :
        adjointFunc(adjointFunc),
        maxActiveVariables(maxActiveVariables),
        maxConstantVariables(maxConstantVariables) {}
  };


  /**
   * @brief A static store for an expression handle.
   *
   * The expression handle is generated from the expression object and stored in a static
   * member variable. Therefore only a pointer needs to be stored on the tape and not the
   * whole expression object.
   *
   * @tparam GradientData  The type of the gradient values that are used in the tape.
   * @tparam         Real  The floating point type of the values that are used in the tape.
   * @tparam    IndexType  The types for the management of the data.
   * @tparam         Expr  The type of the expression from which the handle is generated.
   */
  template<typename GradientData, typename Real, typename IndexType, typename Expr>
  class ExpressionHandleStore {
    private:
      /**
       * @brief The space for the handle.
       */
      static const ExpressionHandle<GradientData, Real, IndexType> handle;
    public:

      /**
       * @brief Retrieve the stored expression handle object.
       *
       * @return The handle for the expression.
       */
      static const ExpressionHandle<GradientData, Real, IndexType>* getHandle() {
        return &handle;
      }
  };

  /**
   * @brief Instantiation of the expression handle store object.
   *
   * @tparam GradientData  The type of the gradient values that are used in the tape.
   * @tparam         Real  The floating point type of the values that are used in the tape.
   * @tparam    IndexType  The types for the management of the data.
   * @tparam         Expr  The type of the expression from which the handle is generated.
   */
  template<typename GradientData, typename Real, typename IndexType, typename Expr>
  const ExpressionHandle<GradientData, Real, IndexType> ExpressionHandleStore<GradientData, Real, IndexType, Expr>::handle(Expr::template evalAdjoint<IndexType, 0, 0>, ExpressionTraits<Expr>::maxActiveVariables, ExpressionTraits<Expr>::maxConstantVariables);
}
