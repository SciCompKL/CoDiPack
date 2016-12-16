
/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015 Chair for Scientific Computing (SciComp), TU Kaiserslautern
 * Homepage: http://www.scicomp.uni-kl.de
 * Contact:  Prof. Nicolas R. Gauger (codi@scicomp.uni-kl.de)
 *
 * Lead developers: Max Sagebaum, Tim Albring (SciComp, TU Kaiserslautern)
 *
 * This file is part of CoDiPack (http://www.scicomp.uni-kl.de/software/codi).
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
 * Authors: Max Sagebaum, Tim Albring, (SciComp, TU Kaiserslautern)
 */

#pragma once

/**
 * @brief Global namespace for CoDiPack - Code Differentiation Package
 */
namespace codi {

  template<typename Real>
  struct InputExpr {

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

  template<typename Real>
  struct CopyExpr {

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
     */
    template<typename IndexType, size_t offset, size_t constantOffset>
    static CODI_INLINE void evalAdjoint(const Real& seed, const IndexType* indices, const PassiveReal* constantValues, const Real* primalValues, Real* adjointValues) {
      CODI_UNUSED(constantValues);
      CODI_UNUSED(primalValues);
      adjointValues[indices[0]] += seed;
    }
  };

  /*
   * @tparam size  The number of arguments of the preaccumulated function.
   */
  template<typename Real, int size>
  struct PreaccExpr {

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

  template<typename Real>
  struct ExpressionTraits<InputExpr<Real> >  {
    static const size_t maxActiveVariables = 0;
    static const size_t maxConstantVariables = 0;
  };

  template<typename Real>
  struct ExpressionTraits<CopyExpr<Real> >  {
    static const size_t maxActiveVariables = 1;
    static const size_t maxConstantVariables = 0;
  };

  template<typename Real, int size>
  struct ExpressionTraits<PreaccExpr<Real, size> >  {
    static const size_t maxActiveVariables = 1;
    static const size_t maxConstantVariables = 0;
  };
}
