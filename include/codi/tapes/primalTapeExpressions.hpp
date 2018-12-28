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

#include "../configure.h"

/**
 * @brief Global namespace for CoDiPack - Code Differentiation Package
 */
namespace codi {

  /**
   * @brief The reverse interpretation of an input operation.
   *
   * @tparam Real  A calculation type that supports all mathematical operations.
   */
  template<typename Real>
  struct InputExpr {

    /** @brief The passive value of the Real type */
    typedef typename TypeTraits<Real>::PassiveReal PassiveReal;

    /**
     * @brief Not implemented for this expression.
     *
     * @param[in]        indices  Not used.
     * @param[in] constantValues  Not used.
     * @param[in]   primalValues  Not used.
     * @return 0.0
     *
     * @tparam          Index  The type for the indices.
     * @tparam         offset  The offset in the index array for the corresponding value.
     * @tparam constantOffset  The offset for the constant values array
     */
    template<typename Index, size_t offset, size_t constantOffset>
    static CODI_INLINE Real getValue(const Index* indices, const PassiveReal* constantValues, const Real* primalValues) {
      CODI_UNUSED(constantValues);

      return primalValues[indices[offset]];
    }

    /**
     * @brief Handle for the input function.
     *
     * It does not need to do anything.
     *
     * @param[in]               seed  The seed for the adjoint of the lhs value.
     * @param[in]            indices  The indices for the arguments of the rhs.
     * @param[in]     constantValues  The array of the constant values in the rhs.
     * @param[in]       primalValues  The global vector with the primal values.
     * @param[in,out]  adjointValues  The global vector with the adjoint values.
     *
     * @tparam          Index  The type for the indices.
     * @tparam  GradientValue  A type that supports add and scalar multiplication.
     * @tparam         offset  The offset in the index array for the corresponding value.
     * @tparam constantOffset  The offset for the constant values array
     */
    template<typename Index, typename GradientValue, size_t offset, size_t constantOffset>
    static CODI_INLINE void evalAdjoint(const PRIMAL_SEED_TYPE& seed, const Index* indices, const PassiveReal* constantValues, const Real* primalValues, PRIMAL_ADJOINT_TYPE* adjointValues) {
      CODI_UNUSED(seed);
      CODI_UNUSED(indices);
      CODI_UNUSED(constantValues);
      CODI_UNUSED(primalValues);
      CODI_UNUSED(adjointValues);
    }

    /**
     * @brief Handle for the input function.
     *
     * Returns the primal value from the input.
     *
     * @param[in]               seed  The seed for the adjoint of the expression.
     * @param[in,out]     lhsAdjoint  The tangent value for the lhs value.
     * @param[in]            indices  The indices for the arguments of the rhs.
     * @param[in]     constantValues  The array of the constant values in the rhs.
     * @param[in]       primalValues  The global vector with the primal values.
     * @param[in,out]  adjointValues  The global vector with the adjoint values.
     *
     * @return The primal value from the input.
     *
     * @tparam          Index  The type for the indices.
     * @tparam  GradientValue  A type that supports add and scalar multiplication.
     * @tparam         offset  The offset in the index array for the corresponding value.
     * @tparam constantOffset  The offset for the constant values array
     */
    template<typename Index, typename GradientValue, size_t offset, size_t constantOffset>
    static CODI_INLINE Real evalTangent(const Real& seed, GradientValue& lhsAdjoint, const Index* indices, const PassiveReal* constantValues, const Real* primalValues, PRIMAL_ADJOINT_TYPE* adjointValues) {
      CODI_UNUSED(seed);
      CODI_UNUSED(lhsAdjoint);
      CODI_UNUSED(constantValues);
      CODI_UNUSED(adjointValues);

      return primalValues[indices[offset]];
    }
  };

  /**
   * @brief The reverse interpretation of a copy operation.
   *
   * @tparam Real  A calculation type that supports all mathematical operations.
   */
  template<typename Real>
  struct CopyExpr {

    /** @brief The passive value of the Real type */
    typedef typename TypeTraits<Real>::PassiveReal PassiveReal;

    /**
     * @brief Returns the value from the primal value vector
     *
     * @param[in]        indices  The index of the rhs
     * @param[in] constantValues  Not used.
     * @param[in]   primalValues  The vector with the primal values.
     * @return 0.0
     *
     * @tparam          Index  The type for the indices.
     * @tparam         offset  The offset in the index array for the corresponding value.
     * @tparam constantOffset  The offset for the constant values array
     */
    template<typename Index, size_t offset, size_t constantOffset>
    static CODI_INLINE Real getValue(const Index* indices, const PassiveReal* constantValues, const Real* primalValues) {
      CODI_UNUSED(constantValues);

      return primalValues[indices[offset]];
    }

    /**
     * @brief Handle for the copy function.
     *
     * It updates the adjoint of the rhs with the seed.
     *
     * @param[in]               seed  The seed for the adjoint of the lhs value.
     * @param[in]            indices  The indices for the arguments of the rhs.
     * @param[in]     constantValues  The array of the constant values in the rhs.
     * @param[in]       primalValues  The global vector with the primal values.
     * @param[in,out]  adjointValues  The global vector with the adjoint values.
     *
     * @tparam          Index  The type for the indices.
     * @tparam  GradientValue  A type that supports add and scalar multiplication.
     * @tparam         offset  The offset in the index array for the corresponding value.
     * @tparam constantOffset  The offset for the constant values array
     */
    template<typename Index, typename GradientValue, size_t offset, size_t constantOffset>
    static CODI_INLINE void evalAdjoint(const PRIMAL_SEED_TYPE& seed, const Index* indices, const PassiveReal* constantValues, const Real* primalValues, PRIMAL_ADJOINT_TYPE* adjointValues) {
      CODI_UNUSED(constantValues);
      CODI_UNUSED(primalValues);
#if CODI_EnableVariableAdjointInterfaceInPrimalTapes
        adjointValues->updateJacobiAdjoint(indices[0], seed);
#else
        adjointValues[indices[0]] += seed;
#endif
    }

    /**
     * @brief Handle for the copy function.
     *
     * Returns the primal value from the input and copies the tangent value.
     *
     * @param[in]               seed  The seed for the adjoint of the expression.
     * @param[in,out]     lhsAdjoint  The tangent value for the lhs value.
     * @param[in]            indices  The indices for the arguments of the rhs.
     * @param[in]     constantValues  The array of the constant values in the rhs.
     * @param[in]       primalValues  The global vector with the primal values.
     * @param[in,out]  adjointValues  The global vector with the adjoint values.
     *
     * @return The primal value from the input.
     *
     * @tparam          Index  The type for the indices.
     * @tparam  GradientValue  A type that supports add and scalar multiplication.
     * @tparam         offset  The offset in the index array for the corresponding value.
     * @tparam constantOffset  The offset for the constant values array
     */
    template<typename Index, typename GradientValue, size_t offset, size_t constantOffset>
    static CODI_INLINE Real evalTangent(const Real& seed, GradientValue& lhsAdjoint, const Index* indices, const PassiveReal* constantValues, const Real* primalValues, PRIMAL_ADJOINT_TYPE* adjointValues) {
      CODI_UNUSED(lhsAdjoint);
      CODI_UNUSED(constantValues);
      CODI_UNUSED(primalValues);

#if CODI_EnableVariableAdjointInterfaceInPrimalTapes
        adjointValues->updateJacobiTangent(indices[offset], seed);
#else
        lhsAdjoint += adjointValues[indices[offset]] * seed;
#endif

      return primalValues[indices[offset]];
    }
  };

  /**
   * @brief The reverse interpretation of a jacobi evaluation.
   *
   * The jacobi values are stored in the passive value vector.
   *
   * @tparam Real  A calculation type that supports all mathematical operations.
   * @tparam size  The number of arguments of the preaccumulated function.
   */
  template<typename Real, size_t size>
  struct PreaccExpr {

    /** @brief The passive value of the Real type */
    typedef typename TypeTraits<Real>::PassiveReal PassiveReal;

    /**
     * @brief Not implemented for this expression.
     *
     * @param[in]        indices  Not used.
     * @param[in] constantValues  Not used.
     * @param[in]   primalValues  Not used.
     * @return 0.0
     *
     * @tparam          Index  The type for the indices.
     * @tparam         offset  The offset in the index array for the corresponding value.
     * @tparam constantOffset  The offset for the constant values array
     */
    template<typename Index, size_t offset, size_t constantOffset>
    static CODI_INLINE Real getValue(const Index* indices, const PassiveReal* constantValues, const Real* primalValues) {
      CODI_UNUSED(indices);
      CODI_UNUSED(constantValues);
      CODI_UNUSED(primalValues);

      std::cerr << "Error: Primal handles are not supported by this expression." << std::endl;
      exit(-1);
      return 0.0;
    }

    /**
     * @brief Handle for the pre accumulation.
     *
     * It assumes that there are size constant values, that represent the Jacobian.
     * Each Jacobian is multiplied with the seed and used to update the rhs value.
     *
     * @param[in]               seed  The seed for the adjoint of the lhs value.
     * @param[in]            indices  The indices for the arguments of the rhs.
     * @param[in]     constantValues  The array of the constant values in the rhs.
     * @param[in]       primalValues  The global vector with the primal values.
     * @param[in,out]  adjointValues  The global vector with the adjoint values.
     *
     * @tparam          Index  The type for the indices.
     * @tparam  GradientValue  A type that supports add and scalar multiplication.
     * @tparam         offset  The offset in the index array for the corresponding value.
     * @tparam constantOffset  The offset for the constant values array
     */
    template<typename Index, typename GradientValue, size_t offset, size_t constantOffset>
    static void evalAdjoint(const PRIMAL_SEED_TYPE& seed, const Index* indices, const PassiveReal* constantValues, const Real* primalValues, PRIMAL_ADJOINT_TYPE* adjointValues) {
      CODI_UNUSED(primalValues);
      CODI_UNUSED(constantValues);

      for(int i = 0; i < (int)size; ++i) {
        // jacobies are stored in the constant values
#if CODI_EnableVariableAdjointInterfaceInPrimalTapes
          adjointValues->updateJacobiAdjoint(indices[i], primalValues[i + 1] * seed);
#else
          adjointValues[indices[i]] += primalValues[i + 1] * seed;
#endif
      }
    }

    /**
     * @brief Currently not supported.
     *
     * @param[in]               seed  The seed for the adjoint of the expression.
     * @param[in,out]     lhsAdjoint  The tangent value for the lhs value.
     * @param[in]            indices  The indices for the arguments of the rhs.
     * @param[in]     constantValues  The array of the constant values in the rhs.
     * @param[in]       primalValues  The global vector with the primal values.
     * @param[in,out]  adjointValues  The global vector with the adjoint values.
     *
     * @return The primal value from the input.
     *
     * @tparam          Index  The type for the indices.
     * @tparam  GradientValue  A type that supports add and scalar multiplication.
     * @tparam         offset  The offset in the index array for the corresponding value.
     * @tparam constantOffset  The offset for the constant values array
     */
    template<typename Index, typename GradientValue, size_t offset, size_t constantOffset>
    static Real evalTangent(const Real& seed, GradientValue& lhsAdjoint, const Index* indices, const PassiveReal* constantValues, const Real* primalValues, PRIMAL_ADJOINT_TYPE* adjointValues) {
      CODI_UNUSED(seed);
      CODI_UNUSED(lhsAdjoint);
      CODI_UNUSED(indices);
      CODI_UNUSED(constantValues);
      CODI_UNUSED(primalValues);
      CODI_UNUSED(adjointValues);

      std::cerr << "Error: Froward handles are not supported by this expression." << std::endl;
      exit(-1);

      return 0.0;
    }
  };

  /**
   * @brief The expression traits for the input expression.
   *
   * @tparam Real  A calculation type that supports all mathematical operations.
   */
  template<typename Real>
  struct ExpressionTraits<InputExpr<Real> >  {
    /** @brief An input expression has no arguments. */
    static const size_t maxActiveVariables = 0;
    /** @brief An input expression has no constant arguments. */
    static const size_t maxConstantVariables = 0;
  };

  /**
   * @brief The expression traits for the copy expression.
   *
   * @tparam Real  A calculation type that supports all mathematical operations.
   */
  template<typename Real>
  struct ExpressionTraits<CopyExpr<Real> >  {
    /** @brief A copy expression has one argument. */
    static const size_t maxActiveVariables = 1;
    /** @brief A copy expression has no constant arguments. */
    static const size_t maxConstantVariables = 0;
  };

  /**
   * @brief The expression traits for the preaccumulation expression.
   *
   * @tparam Real  A calculation type that supports all mathematical operations.
   * @tparam size  The number of arguments of the preaccumulated function.
   */
  template<typename Real, size_t size>
  struct ExpressionTraits<PreaccExpr<Real, size> >  {
    /** @brief The preaccumulation expression has the given size of arguments */
    static const size_t maxActiveVariables = size;
    /** @brief The preaccumulation expression stores the Jacobi entries in the constant value stream.*/
    static const size_t maxConstantVariables = 0;
  };
}
