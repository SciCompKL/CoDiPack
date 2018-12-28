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
 *          Prof. Robin Hogan, (Univ. of Reading).
 *
 * Originally based on Adept 1.0 (http://www.met.rdg.ac.uk/clouds/adept/)
 * released under GPL 3.0 (Copyright (C) 2012-2013 Robin Hogan and the University of Reading).
 */

/*
 * In order to include this file the user has to define the preprocessor macro NAME, FUNCTION and PRIMAL_FUNCTION.
 * NAME contains the name of the generated operation. FUNCTION represents the normal name of that function
 * e.g. 'operator -' or 'sin'. PRIMAL_FUNCTION is the name of a function which can call the primal operator.
 *
 * The defines NAME, FUNCTION and PRIMAL_FUNCTION will be undefined at the end of this template.
 *
 * The user needs to define the derivative computation functions that calculate the derivative
 * with respect to the first and/or second argument. The naming convention is:
 *
 *  dervBB[M]_Name.
 *
 * BB tells which variable is active. Thus we have the combinations
 * 11 -> both are active
 * 10 -> first argument is active
 * 01 -> second argument is active
 *
 * If the M is present the jacobi from the lower expression is not equal to 1.0 and given as an argument.
 * If the M is not present the jacobi is assumed to be equal to 1.0 and not given as an argument.
 *
 * In addition the user needs to define the gradient computations with respect to both arguments.
 * The naming convention is
 *
 *   gradientA_Name
 *   gradientB_Name
 *
 * These functions needs to compute the derivatives with respect to the active variables and
 * call the pushJacobi function on the arguments.That results in the following combinations:
 *
 * derv11_NAME
 * derv11M_NAME
 * derv10_NAME
 * derv10M_NAME
 * derv01_NAME
 * derv01M_NAME
 * gradientA_NAME
 * gradientB_NAME
 */

#ifndef NAME
  #error Please define a name for the binary expression.
#endif
#ifndef FUNCTION
  #error Please define the primal function representation.
#endif
#ifndef PRIMAL_FUNCTION
  #error Please define a function which calls the primal functions representation.
#endif

#include "macros.h"

#define OP NAME
#define OP11 COMBINE(NAME,11)
#define OP10 COMBINE(NAME,10)
#define OP01 COMBINE(NAME,01)
#define FUNC FUNCTION
#define PRIMAL_CALL PRIMAL_FUNCTION
#define DERIVATIVE_FUNC_11   COMBINE(derv11_, NAME)
#define DERIVATIVE_FUNC_11M  COMBINE(derv11M_, NAME)
#define DERIVATIVE_FUNC_10   COMBINE(derv10_, NAME)
#define DERIVATIVE_FUNC_10M  COMBINE(derv10M_, NAME)
#define DERIVATIVE_FUNC_01   COMBINE(derv01_, NAME)
#define DERIVATIVE_FUNC_01M  COMBINE(derv01M_, NAME)

#define GRADIENT_FUNC_A COMBINE(gradientA_, NAME)
#define GRADIENT_FUNC_B COMBINE(gradientB_, NAME)

/* predefine the structs and the functions for higher order derivatives */
template <typename Real, class A, class B> struct OP11;
template <typename Real, class A> struct OP10;
template <typename Real, class B> struct OP01;

template <typename Real, class A, class B>
CODI_INLINE  OP11<Real, A, B> FUNC(const Expression<Real, A>& a, const Expression<Real, B>& b);

template <typename Real, class A>
CODI_INLINE  OP10<Real, A> FUNC(const Expression<Real, A>& a, const typename TypeTraits<Real>::PassiveReal& b);

template <typename Real, class B>
CODI_INLINE  OP01<Real, B> FUNC(const typename TypeTraits<Real>::PassiveReal& a, const Expression<Real, B>& b);

/**
 * @brief Expression implementation for OP with two active variables.
 *
 * @tparam Real The real type used in the active types.
 * @tparam A The expression for the first argument of the function
 * @tparam B The expression for the second argument of the function
 */
template<typename Real, class A, class B>
struct OP11: public Expression<Real, OP11<Real, A, B> > {
  private:

    /** @brief The first argument of the function. */
    CODI_CREATE_STORE_TYPE(A) a_;

    /** @brief The second argument of the function. */
    CODI_CREATE_STORE_TYPE(B) b_;

  public:
    /**
     * @brief The passive type used in the origin.
     *
     * If Real is not an ActiveReal this value corresponds to Real,
     * otherwise the PassiveValue from Real is used.
     */
    typedef typename TypeTraits<Real>::PassiveReal PassiveReal;

    /** @brief Because these are temporary objects they need to be stored as values. */
    static const bool storeAsReference = false;

    /**
     * @brief Stores both arguments of the expression.
     *
     * @param[in] a  First argument of the expression.
     * @param[in] b  Second argument of the expression.
     */
    explicit OP11(const Expression<Real, A>& a, const Expression<Real, B>& b) :
      a_(a.cast()), b_(b.cast()) {}

    /**
     * @brief Calculates the jacobies of the expression and hands them down to the arguments.
     *
     * For f(x,y) it calculates df/dx and df/dy and passes these values as the multipliers to the arguments.
     *
     * @param[in,out] data A helper value which the tape can define and use for the evaluation.
     *
     * @tparam Data The type for the tape data.
     */
    template<typename Data>
    CODI_INLINE void calcGradient(Data& data) const {
#if CODI_DisableCalcGradientSpecialization
        a_.calcGradient(data, GRADIENT_FUNC_A(a_.getValue(), b_.getValue(), getValue()));
        b_.calcGradient(data, GRADIENT_FUNC_B(a_.getValue(), b_.getValue(), getValue()));
#else
        DERIVATIVE_FUNC_11(data, a_, b_, getValue());
#endif
    }

    /**
     * @brief Calculates the jacobies of the expression and hands them down to the arguments.
     *
     * For f(x,y) it calculates multiplier * df/dx and multiplier * df/dy and passes these values as the multipliers to the arguments.
     *
     * @param[in,out]    data A helper value which the tape can define and use for the evaluation.
     * @param[in]  multiplier The Jacobi from the expression where this expression was used as an argument.
     *
     * @tparam Data The type for the tape data.
     */
    template<typename Data>
    CODI_INLINE void calcGradient(Data& data, const Real& multiplier) const {
#if CODI_DisableCalcGradientSpecialization
        a_.calcGradient(data, GRADIENT_FUNC_A(a_.getValue(), b_.getValue(), getValue()) * multiplier);
        b_.calcGradient(data, GRADIENT_FUNC_B(a_.getValue(), b_.getValue(), getValue()) * multiplier);
#else
        DERIVATIVE_FUNC_11M(data, a_, b_, getValue(), multiplier);
#endif
    }

    /**
     * @brief The call is forwarded to both arguments.
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
      a_.pushLazyJacobies(data);
      b_.pushLazyJacobies(data);
    }

    /**
     * @brief Return the numerical value of the expression.
     *
     * @return The value of the expression.
     */
    CODI_INLINE const Real getValue() const {
      return PRIMAL_CALL(a_.getValue(), b_.getValue());
    }

    /**
     * @brief Get the value from a static evaluation context.
     *
     * The method is called in the static evaluation of e.g. a primal value tape.
     * It calls the same method on the arguments with updated offsets for the second argument.
     * The adjustment of the offsets is take from the first argument.
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
    static CODI_INLINE Real getValue(const Index* indices, const PassiveReal* constantValues, const Real* primalValues) {
      const Real aPrimal = A::template getValue<Index, offset, constantOffset>(indices, constantValues, primalValues);
      const Real bPrimal = B::template getValue<Index,
          offset + ExpressionTraits<A>::maxActiveVariables,
          constantOffset + ExpressionTraits<A>::maxConstantVariables>
            (indices, constantValues, primalValues);

      return PRIMAL_CALL(aPrimal, bPrimal);
    }

    /**
     * @brief Calculate the gradient of the expression and update the seed. The updated seed is then
     *        given to the argument expressions.
     *
     * The method is called in the static evaluation of e.g. a primal value tape.
     * It updates the adjoints of the values in the expressions with the calculated
     * adjoint values.
     * It calls the same method on the arguments with updated offsets for the second argument.
     * The adjustment of the offsets is take from the first argument.
     *
     * @param[in]           seed  The seeding for the expression. It is updated in the expressions
     *                            for the operators and used as the update in the terminal points.
     * @param[in]        indices  The indices for the values in the expressions.
     * @param[in] constantValues  The array of constant values in the expression.
     * @param[in]   primalValues  The global primal value vector.
     * @param[in]  adjointValues  The global adjoint value vector.
     *
     * @tparam          Index  The type for the indices.
     * @tparam  GradientValue  A type that supports add and scalar multiplication.
     * @tparam         offset  The offset in the index array for the corresponding value.
     * @tparam constantOffset  The offset for the constant values array
     */
    template<typename Index, typename GradientValue, size_t offset, size_t constantOffset>
    static CODI_INLINE void evalAdjoint(const PRIMAL_SEED_TYPE& seed, const Index* indices,
                                        const PassiveReal* constantValues, const Real* primalValues,
                                        PRIMAL_ADJOINT_TYPE* adjointValues) {
      const Real aPrimal = A::template getValue<Index, offset, constantOffset>(indices, constantValues, primalValues);
      const Real bPrimal = B::template getValue<Index,
          offset + ExpressionTraits<A>::maxActiveVariables,
          constantOffset + ExpressionTraits<A>::maxConstantVariables>
            (indices, constantValues, primalValues);
      const Real resPrimal = PRIMAL_CALL(aPrimal, bPrimal);

      const PRIMAL_SEED_TYPE aJac = GRADIENT_FUNC_A(aPrimal, bPrimal, resPrimal) * seed;
      const PRIMAL_SEED_TYPE bJac = GRADIENT_FUNC_B(aPrimal, bPrimal, resPrimal) * seed;
      A::template evalAdjoint<Index, GradientValue, offset, constantOffset>(aJac, indices, constantValues, primalValues, adjointValues);
      B::template evalAdjoint<Index, GradientValue,
          offset + ExpressionTraits<A>::maxActiveVariables,
          constantOffset + ExpressionTraits<A>::maxConstantVariables>
            (bJac, indices, constantValues, primalValues, adjointValues);
    }

    /**
     * @brief Calculate the gradient of the expression and update the seed. The updated seed is then
     *        given to the argument expressions.
     *
     * The method is called in the static evaluation of e.g. a primal value tape.
     * It computes the tangent direction of the expression with a local reversal of the expression.
     * It calls the same method on the arguments with updated offsets for the second argument.
     * The adjustment of the offsets is take from the first argument.
     *
     * @param[in]           seed  The seeding for the expression. It is updated in the expressions
     *                            for the operators and used as the update in the terminal points.
     * @param[in,out] lhsAdjoint  The tangnet value of the lhs side. This value is updated in the
     *                            arguments of the expression.
     * @param[in]        indices  The indices for the values in the expressions.
     * @param[in] constantValues  The array of constant values in the expression.
     * @param[in]   primalValues  The global primal value vector.
     * @param[in]  adjointValues  The global adjoint value vector.
     *
     * @tparam          Index  The type for the indices.
     * @tparam  GradientValue  A type that supports add and scalar multiplication.
     * @tparam         offset  The offset in the index array for the corresponding value.
     * @tparam constantOffset  The offset for the constant values array
     */
    template<typename Index, typename GradientValue, size_t offset, size_t constantOffset>
    static CODI_INLINE Real evalTangent(const Real& seed, GradientValue& lhsAdjoint, const Index* indices,
                                        const PassiveReal* constantValues, const Real* primalValues,
                                        PRIMAL_ADJOINT_TYPE* adjointValues) {
      const Real aPrimal = A::template getValue<Index, offset, constantOffset>(indices, constantValues, primalValues);
      const Real bPrimal = B::template getValue<Index,
          offset + ExpressionTraits<A>::maxActiveVariables,
          constantOffset + ExpressionTraits<A>::maxConstantVariables>
            (indices, constantValues, primalValues);
      const Real resPrimal = PRIMAL_CALL(aPrimal, bPrimal);

      const Real aJac = GRADIENT_FUNC_A(aPrimal, bPrimal, resPrimal) * seed;
      const Real bJac = GRADIENT_FUNC_B(aPrimal, bPrimal, resPrimal) * seed;
      A::template evalTangent<Index, GradientValue, offset, constantOffset>(aJac, lhsAdjoint, indices, constantValues, primalValues, adjointValues);
      B::template evalTangent<Index, GradientValue,
          offset + ExpressionTraits<A>::maxActiveVariables,
          constantOffset + ExpressionTraits<A>::maxConstantVariables>
            (bJac, lhsAdjoint, indices, constantValues, primalValues, adjointValues);

      return resPrimal;
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
    template<typename Tape, typename Data, typename Func>
    CODI_INLINE void constantValueAction(Tape& tape, Data data, Func func) const {
      a_.constantValueAction(tape, data, func);
      b_.constantValueAction(tape, data, func);
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
      a_.valueAction(data, func);
      b_.valueAction(data, func);
    }
};

/**
 * @brief Expression implementation for OP with one active variables.
 *
 * @tparam Real The real type used in the active types.
 * @tparam A The expression for the first argument of the function
 */
template<typename Real, class A>
struct OP10: public Expression<Real, OP10<Real, A> > {
  private:
    /** @brief The type for the passive values. */
    typedef typename TypeTraits<Real>::PassiveReal PassiveReal;

    /** @brief Active first argument of the function */
    CODI_CREATE_STORE_TYPE(A) a_;

    /** @brief Passive second argument of the function */
    const PassiveReal b_;

  public:

    /** @brief Because these are temporary objects they need to be stored as values. */
    static const bool storeAsReference = false;

    /**
     * @brief Stores both arguments of the expression.
     *
     * @param[in] a  First argument of the expression.
     * @param[in] b  Second argument of the expression.
     */
    explicit OP10(const Expression<Real, A>& a, const PassiveReal& b) :
      a_(a.cast()), b_(b) {}

    /**
     * @brief Calculates the jacobies of the expression and hands them down to the arguments.
     *
     * For f(x,y) it calculates df/dx passes this value as the multiplier to the argument.
     *
     * @param[in,out] data A helper value which the tape can define and use for the evaluation.
     *
     * @tparam Data The type for the tape data.
     */
    template<typename Data>
    CODI_INLINE void calcGradient(Data& data) const {
#if CODI_DisableCalcGradientSpecialization
        a_.calcGradient(data, GRADIENT_FUNC_A(a_.getValue(), b_, getValue()));
#else
        DERIVATIVE_FUNC_10(data, a_, b_, getValue());
#endif
    }

    /**
     * @brief Calculates the jacobies of the expression and hands them down to the arguments.
     *
     * For f(x,y) it calculates multiplier * df/dx and passes this value as the multiplier to the argument.
     *
     * @param[in,out]     data A helper value which the tape can define and use for the evaluation.
     * @param[in]   multiplier The Jacobi from the expression where this expression was used as an argument.
     *
     * @tparam Data The type for the tape data.
     */
    template<typename Data>
    CODI_INLINE void calcGradient(Data& data, const Real& multiplier) const {
#if CODI_DisableCalcGradientSpecialization
        a_.calcGradient(data, GRADIENT_FUNC_A(a_.getValue(), b_, getValue()) * multiplier);
#else
        DERIVATIVE_FUNC_10M(data, a_, b_, getValue(), multiplier);
#endif
    }

    /**
     * @brief The call is forwarded to the active argument.
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
      a_.pushLazyJacobies(data);
    }

    /**
     * @brief Return the numerical value of the expression.
     *
     * @return The value of the expression.
     */
    CODI_INLINE const Real getValue() const {
      return PRIMAL_CALL(a_.getValue(), b_);
    }

    /**
     * @brief Get the value from a static evaluation context.
     *
     * The method is called in the static evaluation of e.g. a primal value tape.
     * It calls the same method on the arguments with updated offsets for the second argument.
     * The adjustment of the offsets is take from the first argument.
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
    static CODI_INLINE Real getValue(const Index* indices, const PassiveReal* constantValues, const Real* primalValues) {
      const Real aPrimal = A::template getValue<Index, offset, constantOffset>(indices, constantValues, primalValues);
      const PassiveReal& bPrimal = constantValues[constantOffset + ExpressionTraits<A>::maxConstantVariables];

      return PRIMAL_CALL(aPrimal, bPrimal);
    }

    /**
     * @brief Calculate the gradient of the expression and update the seed. The updated seed is then
     *        given to the argument expressions.
     *
     * The method is called in the static evaluation of e.g. a primal value tape.
     * It updates the adjoints of the values in the expressions with the calculated
     * adjoint values.
     * It calls the same method on the arguments with updated offsets for the second argument.
     * The adjustment of the offsets is take from the first argument.
     *
     * @param[in]           seed  The seeding for the expression. It is updated in the expressions
     *                            for the operators and used as the update in the terminal points.
     * @param[in]        indices  The indices for the values in the expressions.
     * @param[in] constantValues  The array of constant values in the expression.
     * @param[in]   primalValues  The global primal value vector.
     * @param[in]  adjointValues  The global adjoint value vector.
     *
     * @tparam          Index  The type for the indices.
     * @tparam  GradientValue  A type that supports add and scalar multiplication.
     * @tparam         offset  The offset in the index array for the corresponding value.
     * @tparam constantOffset  The offset for the constant values array
     */
    template<typename Index, typename GradientValue, size_t offset, size_t constantOffset>
    static CODI_INLINE void evalAdjoint(const PRIMAL_SEED_TYPE& seed, const Index* indices, const PassiveReal* constantValues, const Real* primalValues, PRIMAL_ADJOINT_TYPE* adjointValues) {
      const Real aPrimal = A::template getValue<Index, offset, constantOffset>(indices, constantValues, primalValues);
      const PassiveReal& bPrimal = constantValues[constantOffset + ExpressionTraits<A>::maxConstantVariables];
      const Real resPrimal = PRIMAL_CALL(aPrimal, bPrimal);

      const PRIMAL_SEED_TYPE aJac = GRADIENT_FUNC_A(aPrimal, bPrimal, resPrimal) * seed;
      A::template evalAdjoint<Index, GradientValue, offset, constantOffset>(aJac, indices, constantValues, primalValues, adjointValues);
    }

    /**
     * @brief Calculate the gradient of the expression and update the seed. The updated seed is then
     *        given to the argument expressions.
     *
     * The method is called in the static evaluation of e.g. a primal value tape.
     * It computes the tangent direction of the expression with a local reversal of the expression.
     * It calls the same method on the arguments with updated offsets for the second argument.
     * The adjustment of the offsets is take from the first argument.
     *
     * @param[in]           seed  The seeding for the expression. It is updated in the expressions
     *                            for the operators and used as the update in the terminal points.
     * @param[in,out] lhsAdjoint  The tangnet value of the lhs side. This value is updated in the
     *                            arguments of the expression.
     * @param[in]        indices  The indices for the values in the expressions.
     * @param[in] constantValues  The array of constant values in the expression.
     * @param[in]   primalValues  The global primal value vector.
     * @param[in]  adjointValues  The global adjoint value vector.
     *
     * @tparam          Index  The type for the indices.
     * @tparam  GradientValue  A type that supports add and scalar multiplication.
     * @tparam         offset  The offset in the index array for the corresponding value.
     * @tparam constantOffset  The offset for the constant values array
     */
    template<typename Index, typename GradientValue, size_t offset, size_t constantOffset>
    static CODI_INLINE Real evalTangent(const Real& seed, GradientValue& lhsAdjoint, const Index* indices, const PassiveReal* constantValues, const Real* primalValues, PRIMAL_ADJOINT_TYPE* adjointValues) {
      const Real aPrimal = A::template getValue<Index, offset, constantOffset>(indices, constantValues, primalValues);
      const PassiveReal& bPrimal = constantValues[constantOffset + ExpressionTraits<A>::maxConstantVariables];
      const Real resPrimal = PRIMAL_CALL(aPrimal, bPrimal);

      const Real aJac = GRADIENT_FUNC_A(aPrimal, bPrimal, resPrimal) * seed;
      A::template evalTangent<Index, GradientValue, offset, constantOffset>(aJac, lhsAdjoint, indices, constantValues, primalValues, adjointValues);

      return resPrimal;
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
    template<typename Tape, typename Data, typename Func>
    CODI_INLINE void constantValueAction(Tape& tape, Data data, Func func) const {
      a_.constantValueAction(tape, data, func);
      CODI_CALL_MEMBER_FN(tape, func)(data, b_);
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
      a_.valueAction(data, func);
    }
};

/**
 * @brief Expression implementation for OP with one active variables.
 *
 * @tparam Real The real type used in the active types.
 * @tparam B The expression for the second argument of the function
 */
template<typename Real, class B>
struct OP01 : public Expression<Real, OP01<Real, B> > {
  private:

    /** @brief The type for the passive values. */
    typedef typename TypeTraits<Real>::PassiveReal PassiveReal;

    /** @brief Passive first argument of the function */
    const PassiveReal a_;

    /** @brief Active second argument of the function */
    CODI_CREATE_STORE_TYPE(B) b_;
  public:

    /** @brief Because these are temporary objects they need to be stored as values. */
    static const bool storeAsReference = false;

    /**
     * @brief Stores both arguments of the expression.
     *
     * @param[in] a  First argument of the expression.
     * @param[in] b  Second argument of the expression.
     */
    explicit OP01(const PassiveReal& a, const Expression<Real, B>& b) :
      a_(a), b_(b.cast()) {}

    /**
     * @brief Calculates the jacobies of the expression and hands them down to the arguments.
     *
     * For f(x,y) it calculates df/dx passes this value as the multiplier to the argument.
     *
     * @param[in,out] data A helper value which the tape can define and use for the evaluation.
     *
     * @tparam Data The type for the tape data.
     */
    template<typename Data>
    CODI_INLINE void calcGradient(Data& data) const {
#if CODI_DisableCalcGradientSpecialization
        b_.calcGradient(data, GRADIENT_FUNC_B(a_, b_.getValue(), getValue()));
#else
        DERIVATIVE_FUNC_01(data, a_, b_, getValue());
#endif
    }

    /**
     * @brief Calculates the jacobies of the expression and hands them down to the arguments.
     *
     * For f(x,y) it calculates multiplier * df/dx and passes this value as the multiplier to the argument.
     *
     * @param[in,out]    data A helper value which the tape can define and use for the evaluation.
     * @param[in]  multiplier The Jacobi from the expression where this expression was used as an argument.
     *
     * @tparam Data The type for the tape data.
     */
    template<typename Data>
    CODI_INLINE void calcGradient(Data& data, const Real& multiplier) const {
#if CODI_DisableCalcGradientSpecialization
        b_.calcGradient(data, GRADIENT_FUNC_B(a_, b_.getValue(), getValue()) * multiplier);
#else
        DERIVATIVE_FUNC_01M(data, a_, b_, getValue(), multiplier);
#endif
    }

    /**
     * @brief The call is forwarded to the active argument.
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
      b_.pushLazyJacobies(data);
    }

    /**
     * @brief Return the numerical value of the expression.
     *
     * @return The value of the expression.
     */
    CODI_INLINE const Real getValue() const {
      return PRIMAL_CALL(a_, b_.getValue());
    }

    /**
     * @brief Get the value from a static evaluation context.
     *
     * The method is called in the static evaluation of e.g. a primal value tape.
     * It calls the same method on the arguments with updated offsets for the second argument.
     * The adjustment of the offsets is take from the first argument.
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
    static CODI_INLINE Real getValue(const Index* indices, const PassiveReal* constantValues, const Real* primalValues) {
      const PassiveReal& aPrimal = constantValues[constantOffset];
      const Real bPrimal = B::template getValue<Index, offset, constantOffset + 1>(indices, constantValues, primalValues);

      return PRIMAL_CALL(aPrimal, bPrimal);
    }

    /**
     * @brief Calculate the gradient of the expression and update the seed. The updated seed is then
     *        given to the argument expressions.
     *
     * The method is called in the static evaluation of e.g. a primal value tape.
     * It updates the adjoints of the values in the expressions with the calculated
     * adjoint values.
     * It calls the same method on the arguments with updated offsets for the second argument.
     * The adjustment of the offsets is take from the first argument.
     *
     * @param[in]           seed  The seeding for the expression. It is updated in the expressions
     *                            for the operators and used as the update in the terminal points.
     * @param[in]        indices  The indices for the values in the expressions.
     * @param[in] constantValues  The array of constant values in the expression.
     * @param[in]   primalValues  The global primal value vector.
     * @param[in]  adjointValues  The global adjoint value vector.
     *
     * @tparam          Index  The type for the indices.
     * @tparam  GradientValue  A type that supports add and scalar multiplication.
     * @tparam         offset  The offset in the index array for the corresponding value.
     * @tparam constantOffset  The offset for the constant values array
     */
    template<typename Index, typename GradientValue, size_t offset, size_t constantOffset>
    static CODI_INLINE void evalAdjoint(const PRIMAL_SEED_TYPE& seed, const Index* indices, const PassiveReal* constantValues, const Real* primalValues, PRIMAL_ADJOINT_TYPE* adjointValues) {
      const PassiveReal& aPrimal = constantValues[constantOffset];
      const Real bPrimal = B::template getValue<Index, offset, constantOffset + 1>(indices, constantValues, primalValues);
      const Real resPrimal = PRIMAL_CALL(aPrimal, bPrimal);

      const PRIMAL_SEED_TYPE bJac = GRADIENT_FUNC_B(aPrimal, bPrimal, resPrimal) * seed;
      B::template evalAdjoint<Index, GradientValue, offset, constantOffset + 1>(bJac, indices, constantValues, primalValues, adjointValues);
    }

    /**
     * @brief Calculate the gradient of the expression and update the seed. The updated seed is then
     *        given to the argument expressions.
     *
     * The method is called in the static evaluation of e.g. a primal value tape.
     * It computes the tangent direction of the expression with a local reversal of the expression.
     * It calls the same method on the arguments with updated offsets for the second argument.
     * The adjustment of the offsets is take from the first argument.
     *
     * @param[in]           seed  The seeding for the expression. It is updated in the expressions
     *                            for the operators and used as the update in the terminal points.
     * @param[in,out] lhsAdjoint  The tangnet value of the lhs side. This value is updated in the
     *                            arguments of the expression.
     * @param[in]        indices  The indices for the values in the expressions.
     * @param[in] constantValues  The array of constant values in the expression.
     * @param[in]   primalValues  The global primal value vector.
     * @param[in]  adjointValues  The global adjoint value vector.
     *
     * @tparam          Index  The type for the indices.
     * @tparam  GradientValue  A type that supports add and scalar multiplication.
     * @tparam         offset  The offset in the index array for the corresponding value.
     * @tparam constantOffset  The offset for the constant values array
     */
    template<typename Index, typename GradientValue, size_t offset, size_t constantOffset>
    static CODI_INLINE Real evalTangent(const Real& seed, GradientValue& lhsAdjoint, const Index* indices, const PassiveReal* constantValues, const Real* primalValues, PRIMAL_ADJOINT_TYPE* adjointValues) {
      const PassiveReal& aPrimal = constantValues[constantOffset];
      const Real bPrimal = B::template getValue<Index, offset, constantOffset + 1>(indices, constantValues, primalValues);
      const Real resPrimal = PRIMAL_CALL(aPrimal, bPrimal);

      const Real bJac = GRADIENT_FUNC_B(aPrimal, bPrimal, resPrimal) * seed;
      B::template evalTangent<Index, GradientValue, offset, constantOffset + 1>(bJac, lhsAdjoint, indices, constantValues, primalValues, adjointValues);

      return resPrimal;
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
    template<typename Tape, typename Data, typename Func>
    CODI_INLINE void constantValueAction(Tape& tape, Data data, Func func) const {
      CODI_CALL_MEMBER_FN(tape, func)(data, a_);
      b_.constantValueAction(tape, data, func);
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
      b_.valueAction(data, func);
    }
};

/**
 * @brief Specialization of the TypeTraits for the binary operator type.
 *
 * @tparam Real  The floating point value of the active real.
 * @tparam    A  The type of the first argument of the binary operator.
 * @tparam    B  The type of the second argument of the binary operator.
 */
template<typename RealType, typename A, typename B>
class TypeTraits< OP11<RealType, A, B> > {
  public:
    /**
     * @brief The passive type is the passive type of Real.
     */
    typedef typename TypeTraits<RealType>::PassiveReal PassiveReal;

    /**
     * @brief The definition of the Real type for other classes.
     */
    typedef RealType Real;

    /**
     * @brief Get the primal value of the origin of this type.
     * @param[in] t The value from which the primal is extracted.
     * @return The primal value of the origin of this type..
     */
    static const typename TypeTraits<RealType>::PassiveReal getBaseValue(const OP11<RealType, A, B>& t) {
      return TypeTraits<RealType>::getBaseValue(t.getValue());
    }
};

/**
 * @brief Specialization of the TypeTraits for the binary operator type.
 *
 * @tparam Real  The floating point value of the active real.
 * @tparam    A  The type of the first argument of the binary operator.
 */
template<typename RealType, typename A>
class TypeTraits< OP10<RealType, A> > {
  public:
    /**
     * @brief The passive type is the passive type of Real.
     */
    typedef typename TypeTraits<RealType>::PassiveReal PassiveReal;

    /**
     * @brief The definition of the Real type for other classes.
     */
    typedef RealType Real;

    /**
     * @brief Get the primal value of the origin of this type.
     * @param[in] t The value from which the primal is extracted.
     * @return The primal value of the origin of this type..
     */
    static const typename TypeTraits<RealType>::PassiveReal getBaseValue(const OP10<RealType, A>& t) {
      return TypeTraits<RealType>::getBaseValue(t.getValue());
    }
};

/**
 * @brief Specialization of the TypeTraits for the binary operator type.
 *
 * @tparam Real  The floating point value of the active real.
 * @tparam    B  The type of the second argument of the binary operator.
 */
template<typename RealType, typename B>
class TypeTraits< OP01<RealType, B> > {
  public:
    /**
     * @brief The passive type is the passive type of Real.
     */
    typedef typename TypeTraits<RealType>::PassiveReal PassiveReal;

    /**
     * @brief The real is just the type that is defined as the template argument.
     */
    typedef RealType Real;

    /**
     * @brief Get the primal value of the origin of this type.
     * @param[in] t The value from which the primal is extracted.
     * @return The primal value of the origin of this type..
     */
    static const typename TypeTraits<RealType>::PassiveReal getBaseValue(const OP01<RealType, B>& t) {
      return TypeTraits<RealType>::getBaseValue(t.getValue());
    }
};

/**
 * @brief Overload for FUNC with the CoDiPack expressions.
 *
 * @param[in] a  The first argument of the operation.
 * @param[in] b  The second argument of the operation.
 *
 * @return The implementing expression OP.
 *
 * @tparam Real  The real type used in the active types.
 * @tparam    A  The expression for the first argument of the function
 * @tparam    B  The expression for the second argument of the function
 */
template <typename Real, class A, class B>
CODI_INLINE OP11<Real, A, B> FUNC(const Expression<Real, A>& a, const Expression<Real, B>& b) {
  return OP11<Real, A, B>(a.cast(), b.cast());
}
/**
 * @brief Overload for FUNC with the CoDiPack expressions.
 *
 * @param[in] a  The first argument of the operation.
 * @param[in] b  The second argument of the operation.
 *
 * @return The implementing expression OP.
 *
 * @tparam Real  The real type used in the active types.
 * @tparam    A  The expression for the first argument of the function
 */
template <typename Real, class A>
CODI_INLINE OP10<Real, A> FUNC(const Expression<Real, A>& a, const typename TypeTraits<Real>::PassiveReal& b) {
  return OP10<Real, A>(a.cast(), b);
}
/**
 * @brief Overload for FUNC with the CoDiPack expressions.
 *
 * @param[in] a  The first argument of the operation.
 * @param[in] b  The second argument of the operation.
 *
 * @return The implementing expression OP.
 *
 * @tparam Real  The real type used in the active types.
 * @tparam    B  The expression for the second argument of the function
 */
template <typename Real, class B>
CODI_INLINE OP01<Real, B> FUNC(const typename TypeTraits<Real>::PassiveReal& a, const Expression<Real, B>& b) {
  return OP01<Real, B>(a, b.cast());
}

#undef OP
#undef OP11
#undef OP01
#undef OP10
#undef FUNC
#undef PRIMAL_CALL
#undef DERIVATIVE_FUNC_11
#undef DERIVATIVE_FUNC_11M
#undef DERIVATIVE_FUNC_10
#undef DERIVATIVE_FUNC_10M
#undef DERIVATIVE_FUNC_01
#undef DERIVATIVE_FUNC_01M
#undef GRADIENT_FUNC_A
#undef GRADIENT_FUNC_B

#undef PRIMAL_FUNCTION
#undef FUNCTION
#undef NAME
