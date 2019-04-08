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

#include "../activeReal.hpp"
#include "../typeFunctions.hpp"
#include "tapeInterface.hpp"

/**
 * @brief Global namespace for CoDiPack - Code Differentiation Package
 */
namespace codi {
  /**
   * @brief Tape for the tangent or forward AD mode
   *
   * This tape implements the forward or tangent AD mode. For each statement
   * \f[ y = f(x) \f]
   * the rhs of the equation
   * \f[ \dot{y} = \frac{df}{dx}(x)\cdot \dot {x} \f]
   * is evaluated and stored into the GradientData of \f$y\f$. This is done by calling the store routine of the tape
   * in the assignment operator of ActiveReal. Using expression templates, the evaluation of each expression on the rhs leads to an
   * ActiveReal that then calls the pushJacobi routine to
   * add the Jacobian (the partial derivative of the expression with respect to the inputs)
   * multiplied by the tangent value of the input to the tangent value of \f$y\f$.
   *
   * GradientData is just the same as the active type
   * uses for the storage of the floating point values.
   *
   * @tparam           Real  The floating point type of the primal data.
   * @tparam  GradientValue  The floating point type of the tangent data.
   */
  template<typename RealType, typename GradientValueType = RealType>
  class ForwardEvaluation final : public TapeInterface<RealType, GradientValueType, GradientValueType>{
  public:

    /**
     * @brief The real type for the primal values.
     */
    typedef RealType Real;
    /**
     * @brief The real type for the tangent values.
     */
    typedef GradientValueType GradientValue;

    /**
     * @brief The tangent value for the active variable.
     *
     * The tangent data has the same type as the GradientValue data.
     */
    typedef GradientValue GradientData;

    /** @brief Enables code path in CoDiPack that are optimized for Jacobi taping */
    static const bool AllowJacobiOptimization = true;

    /**
     * @brief Evaluates the primal expression and the tangent
     *
     * The store method evaluates the forward AD equation and the primal equation.
     *
     * @param[out]      value  The value of the rhs.
     * @param[out] lhsTangent  The tangent of the lhs.
     * @param[in]         rhs  The expression of the rhs.
     */
    template<typename Rhs>
    CODI_INLINE void store(Real& value, GradientData& lhsTangent, const Rhs& rhs) {
      GradientValue gradient = GradientValue();
      rhs.template calcGradient<GradientValue>(gradient);
      rhs.template pushLazyJacobies<GradientValue>(gradient);
      lhsTangent  = gradient;
      value = rhs.getValue();

#if CODI_AdjointHandle_Tangent
      handleTangentOperation(value, lhsTangent);
#endif

    }

    /**
     * @brief Evaluates the primal expression and the tangent
     *
     * The store method evaluates the forward AD equation and the primal equation.
     *
     * @param[out]      value  The value of the rhs.
     * @param[out] lhsTangent  The tangent of the lhs.
     * @param[in]         rhs  The expression of the rhs.
     */
    CODI_INLINE void store(Real& value, GradientData& lhsTangent, const ActiveReal<ForwardEvaluation<Real> >& rhs) {
      lhsTangent = rhs.getGradient();
      value = rhs.getValue();

#if CODI_AdjointHandle_Tangent
      handleTangentOperation(value, lhsTangent);
#endif
    }

    /**
     * @brief Specialization for store which has a constant value on the rhs.
     *
     * This implementation of store sets the gradient of th active type to zero as the rhs
     * is inactive.
     */
    CODI_INLINE void store(Real& value, GradientData& tangent, const typename TypeTraits<Real>::PassiveReal& rhs) {
      tangent = GradientValue();
      value = rhs;
    }

    /**
     * @brief Adds the jacobi to the tangent value of the expression.
     *
     * This method is called for each value on the rhs. The tangent of the value is added to the
     * tangent of the lhs.
     *
     * @param[in,out] lhsTangent  The tangent of the lhs.
     * @param[in]          value  Not used
     * @param[in]     curTangent  The tangent of the current rhs value.
     *
     * @tparam Data  A Real.
     */
    template<typename Data>
    CODI_INLINE void pushJacobi(Data& lhsTangent, const Real& value, const GradientData& curTangent) {
      CODI_UNUSED(value);
      lhsTangent += curTangent;
    }

    /**
     * @brief Adds the jacobi to the tangent value of the expression.
     *
     * This method is called for each value on the rhs. The tangent of the value times the jacobi is added to the
     * tangent of the lhs.
     *
     * @param[in,out] lhsTangent  The tangent of the lhs.
     * @param[in]         jacobi  The jacobi value of the operation.
     * @param[in]          value  Not used
     * @param[in]     curTangent  The tangent of the current rhs value.
     *
     * @tparam Data  A Real.
     */
    template<typename Data>
    CODI_INLINE void pushJacobi(Data& lhsTangent, const Real& jacobi, const Real& value, const GradientData& curTangent) {
      CODI_UNUSED(value);
      ENABLE_CHECK(OptIgnoreInvalidJacobies, codi::isfinite(jacobi)) {
        lhsTangent += jacobi * curTangent;
      }
    }

    /**
     * @brief Tangent is set to zero.
     *
     * The tangent is initialized with zero.
     *
     * @param[in]    value  Not used.
     * @param[out] tangent  Set to zero.
     */
    CODI_INLINE void initGradientData(Real& value, GradientData& tangent) {
      CODI_UNUSED(value);
      tangent = GradientData();
    }

    /**
     * @brief Nothing to do.
     */
    CODI_INLINE void destroyGradientData(Real& value, GradientData& tangent) {
      CODI_UNUSED(value);
      CODI_UNUSED(tangent);
      /* do nothing */
    }

    /**
     * @brief Calls isTotalZero on the tangent direction.
     *
     * @param[in] gradientData  The tangent direction in this case.
     * @return true if all entries are zero.
     */
    CODI_INLINE bool isGradientTotalZero(const GradientData& gradientData) {
      return codi::isTotalZero(gradientData);
    }

    /**
     * Sets the gradient data to the tangent value.
     *
     * @param[out]   tangent  The tangent value of the active type.
     * @param[in] newTangent  The new tangent value.
     */
    CODI_INLINE void setGradient(GradientData& tangent, const GradientValue& newTangent) {
      tangent = newTangent;
    }

    /**
     * Returns the tangent value of the active type.
     *
     * @param[in]  tangent  The gradient data of the active type is the tangent.
     *
     * @return The tangent value of the active type.
     */
    CODI_INLINE GradientValue getGradient(const GradientData& tangent) const {
      return tangent;
    }

    /**
     * Returns the tangent value of the active type.
     *
     * @param[in]  tangent  The gradient data of the active type is the tangent.
     *
     * @return The reference of the tangent value of the active type.
     */
    CODI_INLINE GradientValue& gradient(GradientData& tangent) {
      return tangent;
    }

    /**
     * Returns the constant tangent value of the active type.
     *
     * @param[in]  tangent  The gradient data of the active type is the tangent.
     *
     * @return The constant reference to the tangent value of the active type.
     */
    CODI_INLINE const GradientValue& gradient(const GradientData& tangent) const {
      return tangent;
    }

    /**
     * @brief Check whether the gradient data is zero.
     *
     * @param[in] tangent  The tangent value for the check.
     * @return False if the gradient data is zero, otherwise returns true.
     */
    bool isActive(const GradientData& tangent) const {
      return !codi::isTotalZero(tangent);
    }
  };
}


