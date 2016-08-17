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
 *          Prof. Robin Hogan, (Univ. of Reading).
 *
 * Originally based on Adept 1.0 (http://www.met.rdg.ac.uk/clouds/adept/)
 * released under GPL 3.0 (Copyright (C) 2012-2013 Robin Hogan and the University of Reading).
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
#ifndef ARG_TYPE
  #error Please define the argument type of the extended argument.
#endif


#include "macros.h"

#define OP NAME
#define FUNC FUNCTION
#define PRIMAL_CALL PRIMAL_FUNCTION
#define GRADIENT_FUNC   COMBINE(grad, NAME)

/* predefine the struct and the function for higher order derivatives */
template<typename Real, class A> struct OP;

template <typename Real, class A>
CODI_INLINE OP<Real, A> FUNC(const Expression<Real, A>& a, ARG_TYPE b);

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
    ARG_TYPE b_;

    /** @brief The result of the function. It is always precomputed. */
    Real result_;
  public:

    /** @brief Because these are temporary objects they need to be stored as values. */
    static const bool storeAsReference = false;

    /**
     * @brief Stores the argument of the expression.
     *
     * @param[in] a Argument of the expression.
     */
    explicit OP(const Expression<Real, A>& a, ARG_TYPE b) :
      a_(a.cast()),
      b_(b),
      result_(PRIMAL_CALL(a.getValue(), b)) {}

  /**
   * @brief Calculates the jacobie of the expression and hands them down to the argument.
   *
   * For f(x) it calculates df/dx and passes this value as the multiplier to the argument.
   *
   * @param[inout] data A helper value which the tape can define and use for the evaluation.
   *
   * @tparam Data The type for the tape data.
   */
  template<typename Data>
  CODI_INLINE void calcGradient(Data& data) const {
    a_.calcGradient(data, GRADIENT_FUNC(a_.getValue(), b_, result_));
  }

  /**
   * @brief Calculates the jacobi of the expression and hands them down to the argument.
   *
   * For f(x) it calculates multiplier * df/dx and passes this value as the multiplier to the argument.
   *
   * @param[inout]     data A helper value which the tape can define and use for the evaluation.
   * @param[in]  multiplier The Jacobi from the expression where this expression was used as an argument.
   *
   * @tparam Data The type for the tape data.
   */
  template<typename Data>
  CODI_INLINE void calcGradient(Data& data, const Real& multiplier) const {
    a_.calcGradient(data, GRADIENT_FUNC(a_.getValue(), b_, result_)*multiplier);
  }

  /**
   * @brief The call is forwarded to the arguments.
   *
   * The method is called for types that accumulate the jacobies before
   * they are pushed to the tape.
   *
   * @param[inout]     data A helper value which the tape can define and use for the evaluation.
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
CODI_INLINE OP<Real, A> FUNC(const Expression<Real, A>& a, ARG_TYPE b) {
  return OP<Real, A>(a.cast(), b);
}

#undef OP
#undef FUNC
#undef PRIMAL_CALL
#undef GRADIENT_FUNC

#undef ARG_TYPE
#undef PRIMAL_FUNCTION
#undef FUNCTION
#undef NAME
