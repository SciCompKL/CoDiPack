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

#include "../../configure.h"
#include "../../evaluateDefinitions.hpp"
#include "../../expressionHandle.hpp"
#include "../../tapeTypes.hpp"
#include "../../typeTraits.hpp"

#include "handleFactoryInterface.hpp"

/**
 * @brief Global namespace for CoDiPack - Code Differentiation Package
 */
namespace codi {

  /**
   * @brief  Storage object for the function pointers for the direct evaluation of an expression.
   *
   * @tparam ReverseTapeTypes  The basic type definitions for the tape. Need to define everything from ReverseTapeTypes.
   */
  template<typename ReverseTapeTypes>
  class FunctionHandle {
    public:

      /**
       * @brief The function pointer to the primal evaluation function
       */
      const typename EvaluateDefinitions<ReverseTapeTypes>::PrimalFunc primalFunc;

      /**
       * @brief The function pointer to the reverse evaluation function
       */
      const typename EvaluateDefinitions<ReverseTapeTypes>::AdjointFunc adjointFunc;

      /**
       * @brief The function pointer to the reverse evaluation function
       */
      const typename EvaluateDefinitions<ReverseTapeTypes>::TangentFunc tangentFunc;

      /**
       * @brief Populate the storae object.
       *
       * @param[in]  primalFunc  The function for the primal evaluation.
       * @param[in] adjointFunc  The function for the reverse evaluation.
       * @param[in] tangentFunc  The function for the tangent evaluation.
       *
       * @tparam P  The type for the primal function object.
       * @tparam A  The type for the reverse function object.
       * @tparam T  The type for the tangent function object.
       */
      template<typename P, typename A, typename T>
      FunctionHandle(const P primalFunc, const A adjointFunc, const T tangentFunc) :
        primalFunc(primalFunc),
        adjointFunc(adjointFunc),
        tangentFunc(tangentFunc) {}
  };


  /**
   * The static store for the function definitions.
   *
   * @tparam Tape  The CoDiPack tape that uses this function store.
   * @tparam Expr  The expression for which this store is instantiated.
   */
  template<typename Tape, typename Expr>
  class FunctionStore {
    private:

      /**
       * The handle that is given to the tape.
       */
      static const FunctionHandle<typename Tape::BaseTypes> handle;
    public:

      /**
       * @brief Retrieve the stored expression handle object.
       *
       * @return The handle for the expression.
       */
      static const FunctionHandle<typename Tape::BaseTypes>* getHandle() {
        return &handle;
      }
  };

  /**
   * The instantiation of the static store object.
   *
   * @tparam Tape  The CoDiPack tape that uses this function store.
   * @tparam Expr  The expression for which this object is instantiated.
   */
  template<typename Tape, typename Expr>
  const FunctionHandle<typename Tape::BaseTypes>
    FunctionStore<Tape, Expr>::handle(
      &Tape::template curryEvaluatePrimalHandle<Expr>,
      &Tape::template curryEvaluateHandle<Expr>,
      &Tape::template curryEvaluateForwardHandle<Expr>);

  /**
   * @brief A factory for function handles, that use static objects to store the data for the function call.
   *
   * The static data of the expression is stored in a static object and a pointer to that
   * static object is returned as the handle.
   *
   * @tparam ReverseTapeTypes  The basic type definitions for a reverse type. Needs to define the
   *                           same types as ReverseTapeTypes
   */
  template<typename ReverseTapeTypes>
  struct StaticFunctionHandleFactory
    // final : public HandleFactoryInterface</* handle type */>
  {
      typedef const FunctionHandle<ReverseTapeTypes>* Handle; /**< Handle type definition */

      /**
       * @brief Create the handle for the given tape and the given expression.
       *
       * @return The pointer to the static object
       *
       * @tparam Expr  The expression that performs the evaluation of the reverse AD operations.
       * @tparam Tape  The tape that is performing the reverse AD evaluation.
       */
      template<typename Expr, typename Tape>
      static CODI_INLINE Handle createHandle() {

        return FunctionStore<Tape, Expr>::getHandle();
      }

      /**
       * @brief The evaluation of the primal handle, that was created by this factory.
       *
       * @param[in]     handle  The handle the was generated by this factory and is called with the arguments.
       * @param[in,out]   args  The other arguments for the function.
       *
       * @tparam Tape  The tape that is performing the reverse AD evaluation.
       * @tparam Args  The arguments for the function.
       */
      template<typename Tape, typename ... Args>
      static CODI_INLINE typename Tape::Real callPrimalHandle(Handle handle, Args&& ... args) {

        return handle->primalFunc(std::forward<Args>(args)...);
      }


      /**
       * @brief The evaluation of the handle, that was created by this factory.
       *
       * The data from the static object is read and used to call the function.
       *
       * @param[in]     handle  The handle the was generated by this factory and is called with the arguments.
       * @param[in,out]   args  The other arguments for the function.
       *
       * @tparam Tape  The tape that is performing the reverse AD evaluation.
       * @tparam Args  The arguments for the function.
       */
      template<typename Tape, typename ... Args>
      static CODI_INLINE void callHandle(Handle handle, Args&& ... args) {

        handle->adjointFunc(std::forward<Args>(args)...);
      }

      /**
       * @brief The evaluation of the forward handle, that was created by this factory.
       *
       * @param[in]     handle  The handle the was generated by this factory and is called with the arguments.
       * @param[in,out]   args  The other arguments for the function.
       *
       * @tparam Tape  The tape that is performing the reverse AD evaluation.
       * @tparam Args  The arguments for the function.
       */
      template<typename Tape, typename ... Args>
      static CODI_INLINE typename Tape::Real callForwardHandle(Handle handle, Args&& ... args) {

        return handle->tangentFunc(std::forward<Args>(args)...);
      }
  };
}
