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

#include "evaluateDefinitions.hpp"
#include "tapeTypes.hpp"
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
   * @tparam ReverseTapeTypes  The basic type definitions for a reverse type. Needs to define the
   *                           same types as ReverseTapeTypes
   */
  template<typename ReverseTapeTypes>
  class ExpressionHandle {
    public:

      /**
       * @brief The function pointer to the primal evaluation function
       */
      const typename EvaluateDefinitions<ReverseTapeTypes>::PrimalExprFunc primalFunc;

      /**
       * @brief The function pointer to the reverse evaluation function
       */
      const typename EvaluateDefinitions<ReverseTapeTypes>::AdjointExprFunc adjointFunc;

      /**
       * @brief The function pointer to the reverse evaluation function
       */
      const typename EvaluateDefinitions<ReverseTapeTypes>::TangentExprFunc tangentFunc;

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
       * @param[in]           primalFunc  The function pointer for the primal evaluation.
       * @param[in]          adjointFunc  The function pointer for the adjoint evaluation.
       * @param[in]          tangentFunc  The function pointer for the tangent evaluation.
       * @param[in]   maxActiveVariables  The number of active variables in the statement.
       * @param[in] maxConstantVariables  The number of constant variables in the statement.
       *
       * @tparam  PrimalFunc  Function with the interface defined in EvaluateDefinitions<ReverseTapeTypes>::PrimalExprFunc
       * @tparam AdjointFunc  Function with the interface defined in EvaluateDefinitions<ReverseTapeTypes>::AdjointExprFunc
       * @tparam TangentFunc  Function with the interface defined in EvaluateDefinitions<ReverseTapeTypes>::TangentExprFunc
       */
      template<typename PrimalFunc, typename AdjointFunc, typename TangentFunc>
      ExpressionHandle(const PrimalFunc primalFunc,
                       const AdjointFunc adjointFunc,
                       const TangentFunc tangentFunc,
                       const size_t maxActiveVariables,
                       const size_t maxConstantVariables) :
        primalFunc(primalFunc),
        adjointFunc(adjointFunc),
        tangentFunc(tangentFunc),
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
   * @tparam Tape  The tape that uses the store. Needs to be a tape type.
   * @tparam Expr  The type of the expression from which the handle is generated.
   */
  template<typename Tape, typename Expr>
  class ExpressionStore {
    private:
      /**
       * @brief The space for the handle.
       */
      static const ExpressionHandle<typename Tape::BaseTypes> handle;
    public:

      /**
       * @brief Retrieve the stored expression handle object.
       *
       * @return The handle for the expression.
       */
      static const ExpressionHandle<typename Tape::BaseTypes>* getHandle() {
        return &handle;
      }
  };

  /**
   * @brief Instantiation of the expression handle store object.
   *
   * @tparam Tape  The tape that uses the store. Needs to be a tape type.
   * @tparam Expr  The type of the expression from which the handle is generated.
   */
  template<typename Tape, typename Expr>
  const ExpressionHandle<typename Tape::BaseTypes> ExpressionStore<Tape, Expr>::handle(
      Expr::template getValue<typename Tape::Index, 0, 0>,
      Expr::template evalAdjoint<typename Tape::Index, typename Tape::GradientValue, 0, 0>,
      Expr::template evalTangent<typename Tape::Index, typename Tape::GradientValue, 0, 0>,
      ExpressionTraits<Expr>::maxActiveVariables,
      ExpressionTraits<Expr>::maxConstantVariables);
}
