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

#include "configure.h"
#include "expressions.hpp"
#include "typeTraits.hpp"
#include "typeFunctions.hpp"
#include "expressionTraits.hpp"

#include <iostream>

/**
 * @brief Global namespace for CoDiPack - Code Differentiation Package
 */
namespace codi {

  /**
   * @brief The overloaded type for the derivative computation.
   *
   * The active type is the starting point for the derivative computation. Through the extension
   * from Expression it can be used in statements like c = a + b to generate expression templates which
   * encapsulate the operation of the right hand side of the statement. During the assignment to the lhs
   * the rhs is evaluated and the results are reported to the tape of this active type.
   * The tape can then decide to store some data or to calculate the derivatives.
   *
   * This type can be nested:
   * \code{.cpp}
   *  typedef ActiveReal<double, Tape1> Real1; // Called origin or start.
   *  typedef ActiveReal<Real1, Tape2> Real2;
   *  typedef ActiveReal<Real2, Tape3> Real3; // Referred to as intermediate types.
   *  ...
   *  typedef ActiveReal<Real(N-1), TapeN> RealN; // Termination type or end of chain
   * \endcode
   *
   * This nesting will also be referred to as a chain. The start of the chain will be named as origin.
   * The origin has the property that its real type is no ActiveReal type. This type will be seen as
   * the passive type. The end of the chain will also be called the termination type.
   *
   * As this type is used to calculate derivative data alongside the programmed and intended
   * computation we call the value of the intended computation the primal computation, primal data or
   * primal expression. Everything corresponding to the derivative computation will most likely be called
   * gradient or derivative.
   *
   * For more information on how to use this class please refer to the #RealForward and #RealReverse.
   *
   * @tparam Tape The tape which handles the derivative calculation. This type has to implement the TapeInterface.
   */
  template<typename Tape>
  class ActiveReal : public Expression<typename Tape::Real, ActiveReal<Tape> > {
  public:

    /**
     * @brief The floating point type used for the calculations.
     */
    typedef typename Tape::Real Real;

    /**
     * @brief Defines that the active real's are stored as references in the expression templates.
     */
    static const bool storeAsReference = true;

    /**
     * @brief Static definition of the tape.
     */
    static Tape globalTape;

    /**
     * @brief The tape used for the derivative calculations.
     */
    typedef Tape TapeType;

    /**
     * @brief The passive type used in the origin.
     *
     * If Real is not an ActiveReal this value corresponds to Real,
     * otherwise the PassiveValue from Real is used.
     */
    typedef typename TypeTraits<Real>::PassiveReal PassiveReal;

    /**
     * @brief The gradient data needed by the tape to store information about the derivatives.
     */
    typedef typename Tape::GradientData GradientData;

    /**
     * @brief The gradient value needed by the tape to compute the derivatives.
     */
    typedef typename Tape::GradientValue GradientValue;

  private:
    /**
     * @brief The primal value of this floating point type.
     */
    Real primalValue;

    /**
     * @brief The gradient data needed by the tape.
     *
     * The active real will never modify this data. It is always passed to the
     * tape as a reference.
     */
    GradientData gradientData;

  public:

    /**
     * @brief Constructs the equivalent of zero and initializes the gradient data.
     */
    CODI_INLINE ActiveReal() : primalValue() {
      globalTape.initGradientData(primalValue, gradientData);
    }

    /**
     * @brief Sets the primal value of the origin and initializes the gradient data.
     *
     * Initializes the value of the start of the ActiveType chain with the value.
     *
     * @param[in] value   The primal value of the active type.
     */
    CODI_INLINE ActiveReal(const PassiveReal& value) : primalValue(value) {
      globalTape.initGradientData(this->primalValue, gradientData);
    }

    /**
     * @brief Sets the primal value of this ActiveReal and sets the gradient after it was initialized.
     *
     * @param[in]    value  The primal value for this type.
     * @param[in] gradient  The gradient value for this type.
     */
    CODI_INLINE ActiveReal(const Real& value, const Real& gradient) : primalValue(value) {
      globalTape.initGradientData(this->primalValue, gradientData);
      globalTape.setGradient(gradientData, gradient);
    }

    /**
     * @brief Forwards the evaluation of the expression to the tape.
     *
     * All logic is handled in the tape. The tape is required to set the primal
     * value to the primal value of the expression.
     *
     * @param[in] rhs The expression which is evaluated.
     *
     * @tparam R The type of the expression.
     */
    template<class R>
    CODI_INLINE ActiveReal(const Expression<Real, R>& rhs) {
      globalTape.initGradientData(this->primalValue, gradientData);
      globalTape.store(primalValue, gradientData, rhs.cast());
    }

    /**
     * @brief Copy Constructor. The logic is handled by the tape.
     *
     * The tape is required to set the primal value of v to the primal value
     * of this type.
     *
     * @param[in] v The value to copy.
     */
    CODI_INLINE ActiveReal(const ActiveReal<Tape>& v) {
      globalTape.initGradientData(primalValue, gradientData);

      if(OptDisableAssignOptimization) {
        *this = PassiveReal(1.0) * v;
      } else {
        globalTape.store(primalValue, gradientData, v);
      }
    }

    /**
     * @brief Call the tape to destroy the gradient data.
     */
    CODI_INLINE ~ActiveReal() {
      globalTape.destroyGradientData(primalValue, gradientData);
    }

    /**
     * @brief Called in the expression evaluation to inform the tape about a partial derivative with the value 1.0.
     *
     * @param[in,out] data A helper value which the tape can define and use for the evaluation.
     *
     * @tparam Data The type for the tape data.
     */
    template<typename Data>
    CODI_INLINE void calcGradient(Data& data) const {
      globalTape.pushJacobi(data, primalValue, gradientData);
    }

    /**
     * @brief Called in the expression evaluation to inform the tape about a partial derivative with the value jacobi.
     *
     * @param[in,out]    data A helper value which the tape can define and use for the evaluation.
     * @param[in]      jacobi The Jacobi from the expression where this expression was used as an argument.
     *
     * @tparam Data The type for the tape data.
     */
    template<typename Data>
    CODI_INLINE void calcGradient(Data& data, const Real& jacobi) const {
      globalTape.pushJacobi(data, jacobi, primalValue, gradientData);
    }

    /**
     * @brief Not needed by this type.
     *
     * The method is called for types that accumulate the jacobies before
     * they are pushed to the tape.
     *
     * @param[in,out]     data A helper value which the tape can define and use for the evaluation.
     *
     * @tparam Data The type for the tape data.
     */
    template<typename Data>
    CODI_INLINE void pushLazyJacobies(Data& data) const {
      CODI_UNUSED(data);
    }

    /**
     * @brief constantValueActions are not executed for active reals.
     *
     * The action is only called for constant values in the expression.
     * E.g. values that have the type PassiveReal.
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
      CODI_UNUSED(tape);
      CODI_UNUSED(data);
      CODI_UNUSED(func);
    }

    /**
     * @brief The action is called on the tape for every active real.
     *
     * @param[in,out] data  The data that can be used by the action.
     * @param[in]     func  The function that is called for every active real in the expression.
     *
     * @tparam     Data  The type of the data for the action.
     * @tparam     Func  The type of the function that is called.
     */
    template<typename Data, typename Func>
    CODI_INLINE void valueAction(Data data, Func func) const {
      CODI_CALL_MEMBER_FN(globalTape, func)(data, primalValue, gradientData);
    }

    /**
     * @brief Helper function for the tape to get its information about this type.
     * @return The gradient data from the tape stored in this type.
     */
    CODI_INLINE GradientData& getGradientData() {
      return gradientData;
    }

    /**
     * @brief Helper function for the tape to get its information about this type.
     * @return The gradient data from the tape stored in this type.
     */
    CODI_INLINE const GradientData& getGradientData() const {
      return gradientData;
    }

    /**
     * @brief Get a reference to the actual gradient value of this instance.
     * @return Reference to the gradient value.
     */
    CODI_INLINE GradientValue& gradient() {
      return globalTape.gradient(gradientData);
    }

    /**
     * @brief Get a constant reference to the actual gradient value of this instance.
     * @return Constant reference to the gradient value.
     */
    CODI_INLINE const GradientValue& gradient() const {
      return globalTape.gradient(gradientData);
    }


    /**
     * @brief Get the value of the gradient of this instance.
     * @return The gradient value.
     */
    CODI_INLINE GradientValue getGradient() const {
      return globalTape.getGradient(gradientData);
    }

    /**
     * @brief Set the value of the gradient of this instance.
     * @param gradient  The new gradient value.
     */
    CODI_INLINE void setGradient(const GradientValue& gradient) {
      globalTape.setGradient(gradientData, gradient);
    }

    /**
     * @brief Check whether the variable is active.
     * @return True if the variables is active, otherwise false.
     */
    CODI_INLINE bool isActive() const{
      return globalTape.isActive(this->gradientData);
    }

    /**
     * @brief Get a reference to the primal value of this instance.
     * @return  Reference to the primal value.
     */
    CODI_INLINE Real& value() {
      return primalValue;
    }

    /**
     * @brief Get a constant reference to the primal value of this instance.
     * @return  Constant reference to the primal value.
     */
    CODI_INLINE const Real& value() const {
      return primalValue;
    }

    /**
     * @brief Get the primal value of this instance.
     * @return The primal value.
     */
    CODI_INLINE const Real& getValue() const {
      return primalValue;
    }

    /**
     * @brief Set the primal value of this instance.
     * @param value The new primal value.
     */
    CODI_INLINE void setValue(const Real& value) {
      this->primalValue = value;
    }

    /**
     * @brief Checks if the primal value and the gradient value is zero.
     * @return true if both are zero.
     */
    CODI_INLINE bool isTotalZero() const {
      return codi::isTotalZero(this->primalValue) && globalTape.isGradientTotalZero(this->gradientData);
    }

    /**
     * @brief Assignment for a passive value on the rhs. E.g a = 1.0;
     *
     * The logic is handled by the tape. The tape is required to set the primal value of
     * the rhs to the primal value of this instance.
     *
     * @param[in] rhs The rhs value.
     * @return Reference to this.
     */
    CODI_INLINE ActiveReal<Tape>& operator=(const PassiveReal& rhs){
      globalTape.store(primalValue, gradientData, rhs);
      return *this;
    }

    /**
     * @brief Assignment for an expression on the rhs. E.g a = x + y;
     *
     * The logic is handled by the tape. The tape is required to set the primal value of
     * the rhs to the primal value of this instance.
     *
     * @param[in] rhs The expression on the rhs.
     * @return Reference to this.
     */
    template<class R>
    CODI_INLINE ActiveReal<Tape>& operator=(const Expression<Real, R>& rhs){
      globalTape.store(primalValue, gradientData, rhs.cast());
      return *this;
    }

    /**
     * @brief Assignment for an ActiveReal on the rhs. E.g a = x;
     *
     * The logic is handled by the tape. The tape is required to set the primal value of
     * the rhs to the primal value of this instance.
     *
     * @param[in] rhs The other value on the rhs.
     * @return Reference to this.
     */
    CODI_INLINE ActiveReal<Tape>& operator=(const ActiveReal<Tape>& rhs) {
      if(OptDisableAssignOptimization) {
           *this = PassiveReal(1.0) * rhs;
      } else {
          globalTape.store(primalValue, gradientData, rhs);
      }
      return *this;
    }

    /**
     * @brief The expression is unfolded to *this = *this + rhs
     *
     * @param[in] rhs The expression on the rhs.
     * @tparam R The type of the expression on the rhs.
     */
    template<class R>
    CODI_INLINE ActiveReal<Tape>& operator+=(const Expression<Real, R>& rhs) {
      return *this = (*this + rhs);
    }
    /**
     * @brief The expression is unfolded to *this = *this - rhs
     *
     * @param[in] rhs The expression on the rhs.
     * @tparam R The type of the expression on the rhs.
     */
    template<class R>
    CODI_INLINE ActiveReal<Tape>& operator-=(const Expression<Real, R>& rhs) {
      return *this = (*this - rhs);
    }
    /**
     * @brief The expression is unfolded to *this = *this * rhs
     *
     * @param[in] rhs The expression on the rhs.
     * @tparam R The type of the expression on the rhs.
     */
    template<class R>
    CODI_INLINE ActiveReal<Tape>& operator*=(const Expression<Real, R>& rhs) {
      return *this = (*this * rhs);
    }
    /**
     * @brief The expression is unfolded to *this = *this / rhs
     *
     * @param[in] rhs The expression on the rhs.
     * @tparam R The type of the expression on the rhs.
     */
    template<class R>
    CODI_INLINE ActiveReal<Tape>& operator/=(const Expression<Real, R>& rhs) {
      return *this = (*this / rhs);
    }

    /**
     * @brief Optimization for expressions like a += 3.0;
     *
     * For this expression the derivative value is not modified.
     *
     * This logic is hidden from the tape.
     *
     * @param[in] rhs The passive value on the rhs.
     */
    CODI_INLINE ActiveReal<Tape>& operator+=(const PassiveReal& rhs) {
      if(Tape::AllowJacobiOptimization) {
        primalValue += rhs;
      } else {
        *this = (*this + rhs);
      }
      return *this;
    }
    /**
     * @brief Optimization for expressions like a -= 3.0;
     *
     * For this expression the derivative value is not modified.
     *
     * This logic is hidden from the tape.
     *
     * @param[in] rhs The passive value on the rhs.
     */
    CODI_INLINE ActiveReal<Tape>& operator-=(const PassiveReal& rhs) {
      if(Tape::AllowJacobiOptimization) {
        primalValue -= rhs;
      } else {
        *this = (*this - rhs);
      }
      return *this;
    }
    /**
     * @brief The expression is unfolded to *this = *this * rhs
     *
     * @param[in] rhs The passive value on the rhs.
     */
    CODI_INLINE ActiveReal<Tape>& operator*=(const PassiveReal& rhs) {
      return *this = (*this * rhs);
    }
    /**
     * @brief The expression is unfolded to *this = *this * rhs
     *
     * @param[in] rhs The passive value on the rhs.
     */
    CODI_INLINE ActiveReal<Tape>& operator/=(const PassiveReal& rhs) {
      return *this = (*this / rhs);
    }

    /**
     * @brief The expression is unfolded to *this += 1.0
     */
    CODI_INLINE ActiveReal<Tape> operator++() {
      return *this = *this + PassiveReal(1.0);
    }

    /**
     * @brief The expression is unfolded to *this += 1.0
     *
     * @param u Indicator for postfix operator.
     */
    CODI_INLINE ActiveReal<Tape> operator++(int u) {
      CODI_UNUSED(u);
      ActiveReal<Tape> r(*this);
      *this = *this + PassiveReal(1.0);
      return r;
    }
    /**
     * @brief The expression is unfolded to *this -= 1.0
     */
    CODI_INLINE ActiveReal<Tape> operator--() {
      return *this = *this - PassiveReal(1.0);
    }
    /**
     * @brief The expression is unfolded to *this -= 1.0
     *
     * @param u Indicator for postfix operator.
     */
    CODI_INLINE ActiveReal<Tape> operator--(int u) {
      CODI_UNUSED(u);
      ActiveReal<Tape> r(*this);
      *this = *this - PassiveReal(1.0);
      return r;
    }

    /**
     * @brief Get the reference to the global tape for this type.
     * @return  The global reference to the tape.
     */
    static CODI_INLINE Tape& getGlobalTape() {
      return globalTape;
    }

    /**
     * @brief Get the value from a static evaluation context.
     *
     * The method is called in the static evaluation of e.g. a primal value tape.
     * It returns the primal value the is defined by the corresponding index.
     * The index is defined by the offset in the index vector.
     *
     * E.g. primalValues[indices[offset]]
     *
     * @param[in]        indices  The indices for the values in the expressions.
     * @param[in] constantValues  The array of constant values in the expression.
     * @param[in]   primalValues  The global primal value vector.
     *
     * @return The corresponding primal value for the active real.
     *
     * @tparam          Index  The type for the indices.
     * @tparam         offset  The offset in the index array for the corresponding value.
     * @tparam constantOffset  The offset for the constant values array
     */
    template<typename Index, size_t offset, size_t constantOffset>
    static CODI_INLINE const Real& getValue(const Index* indices, const PassiveReal* constantValues, const Real* primalValues) {
      CODI_UNUSED(constantValues);
      return primalValues[indices[offset]];
    }

    /**
     * @brief Update the adjoint of the corresponding value in the expression.
     *
     * The method is called in the static evaluation of e.g. a primal value tape.
     * It updates the adjoints of the values in the expressions with the calculated
     * adjoint values.
     *
     * @param[in]           seed  The seeding for the expression. It is updated in the expressions
     *                            for the operators and used as the update in the terminal points.
     * @param[in]        indices  The indices for the values in the expressions.
     * @param[in] constantValues  The array of constant values in the expression.
     * @param[in]   primalValues  The global primal value vector.
     * @param[in]  adjointValues  The global adjoint value vector.
     *
     * @tparam          Index  The type for the indices.
     * @tparam  GradientValue  The type for the gradient values. It needs to provide add functions and a scalar copy.
     * @tparam         offset  The offset in the index array for the corresponding value.
     * @tparam constantOffset  The offset for the constant values array
     */
    template<typename Index, typename GradientValue, size_t offset, size_t constantOffset>
    static CODI_INLINE void evalAdjoint(const PRIMAL_SEED_TYPE& seed, const Index* indices,
                                        const PassiveReal* constantValues, const Real* primalValues,
                                        PRIMAL_ADJOINT_TYPE* adjointValues) {
      CODI_UNUSED(constantValues);
      CODI_UNUSED(primalValues);

      ENABLE_CHECK(OptIgnoreInvalidJacobies, codi::isfinite(seed)) {
#if CODI_EnableVariableAdjointInterfaceInPrimalTapes
        adjointValues->updateJacobiAdjoint(indices[offset], seed);
#else
        adjointValues[indices[offset]] += seed;
#endif
      }
    }

    /**
     * @brief Add the tangent influence of this value in the expression.
     *
     * The method is called in the static evaluation of e.g. a primal value tape.
     * It updates the tangent value of the expression with the calculated
     * adjoint value from the reversal of this expression.
     *
     * @param[in]           seed  The seeding for the expression. It is updated in the expressions
     *                            for the operators and used as the update in the terminal points.
     * @param[in,out] lhsAdjoint  The tangent value for the lhs.
     * @param[in]        indices  The indices for the values in the expressions.
     * @param[in] constantValues  The array of constant values in the expression.
     * @param[in]   primalValues  The global primal value vector.
     * @param[in]  adjointValues  The global adjoint value vector.
     *
     * @tparam          Index  The type for the indices.
     * @tparam  GradientValue  The type for the gradient values. It needs to provide add functions and a scalar copy.
     * @tparam         offset  The offset in the index array for the corresponding value.
     * @tparam constantOffset  The offset for the constant values array
     */
    template<typename Index, typename GradientValue, size_t offset, size_t constantOffset>
    static CODI_INLINE Real evalTangent(const Real& seed, GradientValue& lhsAdjoint, const Index* indices,
                                        const PassiveReal* constantValues, const Real* primalValues,
                                        PRIMAL_ADJOINT_TYPE* adjointValues) {
      CODI_UNUSED(lhsAdjoint);
      CODI_UNUSED(constantValues);
      CODI_UNUSED(primalValues);

      ENABLE_CHECK(OptIgnoreInvalidJacobies, codi::isfinite(seed)) {
#if CODI_EnableVariableAdjointInterfaceInPrimalTapes
        adjointValues->updateJacobiTangent(indices[offset], seed);
#else
        lhsAdjoint += adjointValues[indices[offset]] * seed;
#endif
      }

      return primalValues[indices[offset]];
    }
  };

  /**
   * @brief Specialization of the TypeTraits for the ActiveReal type.
   *
   * @tparam Tape The tape of the active real.
   */
  template<typename Tape>
  class TypeTraits<ActiveReal<Tape> > {
    public:

      /** @brief The the calculation type. */
      typedef typename Tape::Real Real;

      /** @brief The passive type is the passive type of Real. */
      typedef typename TypeTraits<Real>::PassiveReal PassiveReal;

      /** @brief The maximum derivative order that the active type contains. */
      static const size_t MaxDerivativeOrder = 1 + TypeTraits<Real>::MaxDerivativeOrder;

      /**
       * @brief Get the primal value of the origin of this type.
       * @param[in] t The value from which the primal is extracted.
       * @return The primal value of the origin of this type..
       */
      static const PassiveReal getBaseValue(const ActiveReal<Tape>& t) {
        return TypeTraits<Real>::getBaseValue(t.getValue());
      }
  };

  /**
   * @brief The instantiation of the tape for the ActiveReal.
   */
  template<typename Tape>
  Tape ActiveReal<Tape>::globalTape;

  /**
   * @brief Specialization of the ExpressionTraits for the ActiveReal type.
   *
   * @tparam Tape The tape of the active real.
   */
  template<typename Tape>
  struct ExpressionTraits<ActiveReal<Tape> >  {
    /**
     * @brief The maximum number of active values for an ActiveReal is one.
     */
    static const size_t maxActiveVariables = 1;

    /**
     * @brief The maximum number of passive values for an ActiveReal is zero.
     */
    static const size_t maxConstantVariables = 0;
  };

  /**
   * @brief The primal value of the origin is written to the stream.
   *
   * @param[in,out] os The stream for the operation.
   * @param[in]    rhs The expression on the rhs of the stream operation.
   *
   * @return The modified stream.
   *
   * @tparam Real The floating point value of the active real.
   * @tparam Tape The tape of the active real.
   */
  template<typename Real, class R>
  std::ostream& operator<<(std::ostream& os, const Expression<Real, R>& rhs){
    os << rhs.getValue();
    return os;
  }

  /**
   * @brief A passive value is read from the stream and set to the primal value of the origin.
   *
   * @param[in,out] os The stream for the operation.
   * @param[out]   rhs The activeValue on the rhs of the stream operation.
   *
   * @return The modified stream.
   *
   * @tparam Tape The tape of the active real.
   */
  template<typename Tape>
  std::istream& operator>>(std::istream& os, ActiveReal<Tape>& rhs){
    typename Tape::Real temp;
    os >> temp;
    rhs.setValue(temp);
    return os;
  }
}
