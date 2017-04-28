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
     * @tparam      IndexType  The type for the indices.
     * @tparam         offset  The offset in the index array for the corresponding value.
     * @tparam constantOffset  The offset for the constant values array
     */
    template<typename IndexType, size_t offset, size_t constantOffset>
    static CODI_INLINE void evalAdjoint(const Real& seed, const IndexType* indices, const PassiveReal* constantValues, const Real* primalValues, Real* adjointValues) {
      CODI_UNUSED(seed);
      CODI_UNUSED(indices);
      CODI_UNUSED(constantValues);
      CODI_UNUSED(primalValues);
      CODI_UNUSED(adjointValues);
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
     * @tparam      IndexType  The type for the indices.
     * @tparam         offset  The offset in the index array for the corresponding value.
     * @tparam constantOffset  The offset for the constant values array
     */
    template<typename IndexType, size_t offset, size_t constantOffset>
    static CODI_INLINE void evalAdjoint(const Real& seed, const IndexType* indices, const PassiveReal* constantValues, const Real* primalValues, Real* adjointValues) {
      CODI_UNUSED(constantValues);
      CODI_UNUSED(primalValues);
      adjointValues[indices[0]] += seed;
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
  template<typename Real, int size>
  struct PreaccExpr {

    /** @brief The passive value of the Real type */
    typedef typename TypeTraits<Real>::PassiveReal PassiveReal;

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
     * @tparam      IndexType  The type for the indices.
     * @tparam         offset  The offset in the index array for the corresponding value.
     * @tparam constantOffset  The offset for the constant values array
     */
    template<typename IndexType, size_t offset, size_t constantOffset>
    static void evalAdjoint(const Real& seed, const IndexType* indices, const PassiveReal* constantValues, const Real* primalValues, Real* adjointValues) {
      CODI_UNUSED(primalValues);
      for(int i = 0; i < size; ++i) {
        // jacobies are stored in the constant values
        adjointValues[indices[i]] += constantValues[i] * seed;
      }
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
    static const size_t maxConstantVariables = size;
  };
}
