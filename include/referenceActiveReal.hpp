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

#include "activeReal.hpp"
#include "expressions.hpp"
#include "typeTraits.hpp"
#include "expressionTraits.hpp"
#include <iostream>

namespace codi {

  /**
   * @brief A helper structure for optimizing the case when a variable is used multiple times in an expression.
   *
   * If in a statement like
   * \code{.cpp}
   *    y = sin(x)*cos(x);
   * \endcode
   * the ReferenceActiveReal is used  like
   * \code{.cpp}
   *    RefReal xRef = x;
   *    y = sin(xRef)*cos(xRef);
   * \endcode
   * then only one argument is stored instead of two.
   *
   * For a more detailed explanation see @ref Tutorial5
   *
   * @tparam ActiveType  The type that is referenced in this structure.
   */
  template<typename ActiveType>
  class ReferenceActiveReal : public Expression<typename ActiveType::Real, ReferenceActiveReal<ActiveType> > {
  public:

    /**
     * @brief This type needs to be stored as a reference.
     */
    static const bool storeAsReference = true;

    /**
     * @brief The tape type for other classes.
     */
    typedef typename ActiveType::TapeType TapeType;

    /**
     * @brief The underlying floating point type for other classes.
     */
    typedef typename ActiveType::Real Real;

    /**
     * @brief The passive floating point type for other classes.
     */
    typedef typename TypeTraits<Real>::PassiveReal PassiveReal;

    /**
     * @brief The value of the gradient data for other classes.
     */
    typedef typename TapeType::GradientData GradientData;

  private:

    /**
     * @brief The reference to the value for which the optimization is done.
     */
    const ActiveType& reference;

    /**
     * @brief The accumulated Jacobian.
     *
     * The member is declared mutable because the arguments the expressions are usually constant.
     */
    mutable Real jacobi;

  public:

    /**
     * @brief Construct a ReferenceActiveReal the accumulates the Jacobies for the referenced ActiveReal.
     *
     * @param[in] reference  The value for which the accumulation is done.
     */
    CODI_INLINE ReferenceActiveReal(const ActiveType& reference) :
      reference(reference),
      jacobi() {}


    /**
     * @brief The call is not forwarded to the tape. Instead the local Jacobian is updated.
     *
     * @param[in,out] data A helper value which the tape can define and use for the evaluation.
     *
     * @tparam Data The type for the tape data.
     */
    template<typename Data>
    CODI_INLINE void calcGradient(Data& data) const {
      this->jacobi += 1.0;
    }

    /**
     * @brief The call is not forwarded to the tape. Instead the local Jacobian is updated.
     *
     * @param[in,out] data  A helper value which the tape can define and use for the evaluation.
     * @param[in]   jacobi  The value of the Jacobian with respect to this argument of the expression.
     *
     * @tparam Data The type for the tape data.
     */
    template<typename Data>
    CODI_INLINE void calcGradient(Data& data, const Real& jacobi) const {
      this->jacobi += jacobi;
    }

    /**
     * @brief The method calcGradient is called on the referenced value.
     *
     * The argument for the Jacobian is the accumulated Jacobian that is stored in this structure.
     *
     * @param[in,out] data A helper value which the tape can define and use for the evaluation.
     *
     * @tparam Data The type for the tape data.
     */
    template<typename Data>
    CODI_INLINE void pushLazyJacobies(Data& data) const {
      if(0.0 != jacobi) {
        reference.calcGradient(data, jacobi);
        jacobi = 0.0; // reset jacobi for the next statement or the next call for this statement
      }
    }

    /**
     * @brief Returns the gradient data from the referenced ActiveReal.
     *
     * @return The gradient data from the referenced ActiveReal.
     */
    CODI_INLINE const GradientData& getGradientData() const {
      return reference.getGradientData();
    }

    /**
     * @brief Returns the gradient from the referenced ActiveReal.
     *
     * @return The gradient from the referenced ActiveReal.
     */
    CODI_INLINE Real getGradient() const {
      return reference.getGradient();
    }

    /**
     * @brief Returns the value from the referenced ActiveReal.
     *
     * @return The value from the referenced ActiveReal.
     */
    CODI_INLINE const Real& getValue() const {
      return reference.getValue();
    }

    /**
     * @brief Get the value from a static evaluation context.
     *
     * The call is forwarded to the referenced ActiveReal type.
     *
     * @param[in]        indices  The indices for the values in the expressions.
     * @param[in] constantValues  The array of constant values in the expression.
     * @param[in]   primalValues  The global primal value vector.
     *
     * @return The corresponding primal value for the active real.
     *
     * @tparam      IndexType  The type for the indices.
     * @tparam         offset  The offset in the index array for the corresponding value.
     * @tparam constantOffset  The offset for the constant values array
     */
    template<typename IndexType, size_t offset, size_t constantOffset>
    static CODI_INLINE const Real& getValue(const IndexType* indices, const PassiveReal* constantValues, const Real* primalValues) {
      return ActiveType::template getValue<IndexType, offset, constantOffset>(indices, constantValues, primalValues);
    }


    /**
     * @brief Calculate the gradient of the expression and update the seed. The updated seed is then
     *        given to the argument expressions.
     *
     * The call is forwarded to the referenced ActiveReal type.
     *
     * @param[in]           seed  The seeding for the expression. It is updated in the expressions
     *                            for the operators and used as the update in the terminal points.
     * @param[in]        indices  The indices for the values in the expressions.
     * @param[in] constantValues  The array of constant values in the expression.
     * @param[in]   primalValues  The global primal value vector.
     * @param[in]  adjointValues  The global adjoint value vector.
     *
     * @tparam      IndexType  The type for the indices.
     * @tparam         offset  The offset in the index array for the corresponding value.
     * @tparam constantOffset  The offset for the constant values array
     */
    template<typename IndexType, size_t offset, size_t constantOffset>
    static CODI_INLINE void evalAdjoint(const Real& seed, const IndexType* indices, const PassiveReal* constantValues, const Real* primalValues, Real* adjointValues) {
      ActiveType::template evalAdjoint<IndexType, offset, constantOffset>(seed, indices, constantValues, primalValues, adjointValues);
    }

    /**
     * @brief The action is forwarded to the referenced ActiveReal.
     *
     * @param[in,out] data  The data that can be used by the action.
     * @param[in]     func  The function that is called for every active real in the expression.
     *
     * @tparam     Data  The type of the data for the action.
     * @tparam     Func  The type of the function that is called.
     */
    template<typename Data, typename Func>
    CODI_INLINE void valueAction(Data data, Func func) const {
      reference.valueAction(data, func);
    }

    /**
     * @brief The action is forwarded to the referenced ActiveReal.
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
      reference.constantValueAction(tape, data, func);
    }

  private:
    /**
     * @brief Forbid assignment onn this type.
     */
    CODI_INLINE ReferenceActiveReal<ActiveType>& operator=(const ReferenceActiveReal& rhs){}
  };

  /**
   * @brief Specialization of the TypeTraits for the ReferenceActiveReal type.
   *
   * @tparam ActiveType  The active type which is stored in this reference object.
   */
  template<typename ActiveType>
  class TypeTraits<ReferenceActiveReal<ActiveType> > {
  public:

    /**
     * @brief The tape type for other classes.
     */
    typedef typename ActiveType::TapeType Tape;

    /**
     * @brief The the calculation type.
     */
    typedef typename Tape::Real Real;

    /**
     * @brief The passive type is the passive type of Real.
     */
    typedef typename TypeTraits<Real>::PassiveReal PassiveReal;

    /**
     * @brief Get the primal value of the origin of this type.
     * @param[in] t The value from which the primal is extracted.
     * @return The primal value of the origin of this type..
     */
    static const typename TypeTraits<Real>::PassiveReal getBaseValue(const ReferenceActiveReal<ActiveType>& t) {
      return TypeTraits<Real>::getBaseValue(t.getValue());
    }
  };

  /**
   * @brief Specialization of the ExpressionTraits for the ReferenceActiveReal type.
   *
   * @tparam ActiveType  The active type which is stored in this reference object.
   */
  template<typename ActiveType>
  struct ExpressionTraits<ReferenceActiveReal<ActiveType> >  {
    /**
     * @brief The maximum number of active values for an ReferenceActiveReal is one.
     */
    static const size_t maxActiveVariables = 1;

    static const size_t maxConstantVariables = 0;
  };
}
