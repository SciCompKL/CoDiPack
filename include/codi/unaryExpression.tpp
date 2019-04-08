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
 * The user needs to define the following functions:
 *
 * gradNAME: Computes the derivative df(x)/dx for y = f(x)
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
#define FUNC FUNCTION
#define PRIMAL_CALL PRIMAL_FUNCTION
#define GRADIENT_FUNC   COMBINE(gradient_, NAME)

/* predefine the struct and the function for higher order derivatives */
template<typename Real, class A> struct OP;

template <typename Real, class A>
CODI_INLINE OP<Real, A> FUNC(const Expression<Real, A>& a);

/**
 * @brief Expression implementation for OP.
 *
 * @tparam Real  The real type used in the active types.
 * @tparam    A  The expression for the argument of the function.
 */
template<typename Real, class A>
struct OP : public Expression<Real, OP<Real, A> > {
  private:

    /** @brief The argument of the function. */
    CODI_CREATE_STORE_TYPE(A) a_;

    /** @brief The result of the function. It is always precomputed. */
    Real result_;
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
     * @brief Stores the argument of the expression.
     *
     * @param[in] a Argument of the expression.
     */
    explicit OP(const Expression<Real, A>& a) :
      a_(a.cast()),
      result_(PRIMAL_CALL(a.getValue())) {}

  /**
   * @brief Calculates the jacobie of the expression and hands them down to the argument.
   *
   * For f(x) it calculates df/dx and passes this value as the multiplier to the argument.
   *
   * @param[in,out] data A helper value which the tape can define and use for the evaluation.
   *
   * @tparam Data The type for the tape data.
   */
  template<typename Data>
  CODI_INLINE void calcGradient(Data& data) const {
    a_.calcGradient(data, GRADIENT_FUNC(a_.getValue(), result_));
  }

  /**
   * @brief Calculates the jacobi of the expression and hands them down to the argument.
   *
   * For f(x) it calculates multiplier * df/dx and passes this value as the multiplier to the argument.
   *
   * @param[in,out]     data A helper value which the tape can define and use for the evaluation.
   * @param[in]   multiplier The Jacobi from the expression where this expression was used as an argument.
   *
   * @tparam Data The type for the tape data.
   */
  template<typename Data>
  CODI_INLINE void calcGradient(Data& data, const Real& multiplier) const {
    a_.calcGradient(data, GRADIENT_FUNC(a_.getValue(), result_)*multiplier);
  }

  /**
   * @brief The call is forwarded to the arguments.
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
  CODI_INLINE const Real& getValue() const {
    return result_;
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

    return PRIMAL_CALL(aPrimal);
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
   *                           for the operators and used as the update in the terminal points.
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
  static CODI_INLINE void evalAdjoint(const PRIMAL_SEED_TYPE& seed, const Index* indices, const PassiveReal* constantValues, const Real* primalValues, PRIMAL_ADJOINT_TYPE* adjointValues) {
    const Real aPrimal = A::template getValue<Index, offset, constantOffset>(indices, constantValues, primalValues);
    const Real resPrimal = PRIMAL_CALL(aPrimal);

    const PRIMAL_SEED_TYPE aJac = GRADIENT_FUNC(aPrimal, resPrimal) * seed;
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
    const Real resPrimal = PRIMAL_CALL(aPrimal);

    const Real aJac = GRADIENT_FUNC(aPrimal, resPrimal) * seed;
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
 * @brief Specialization of the TypeTraits for the unary operator type.
 *
 * @tparam Real  The floating point value of the active real.
 * @tparam    A  The type of the argument of the unary operator.
 */
template<typename RealType, typename A>
class TypeTraits< OP<RealType, A> > {
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
    static const typename TypeTraits<RealType>::PassiveReal getBaseValue(const OP<RealType, A>& t) {
      return TypeTraits<RealType>::getBaseValue(t.getValue());
    }
};

/**
 * @brief Overload for FUNC with the CoDiPack expressions.
 *
 * @param[in] a The argument of the operation.
 *
 * @return The implementing expression OP.
 *
 * @tparam Real The real type used in the active types.
 * @tparam A The expression for the first argument of the function.
 */
template <typename Real, class A>
CODI_INLINE OP<Real, A> FUNC(const Expression<Real, A>& a) {
  return OP<Real, A>(a.cast());
}

#undef OP
#undef FUNC
#undef PRIMAL_CALL
#undef GRADIENT_FUNC

#undef PRIMAL_FUNCTION
#undef FUNCTION
#undef NAME
