/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2020 Chair for Scientific Computing (SciComp), TU Kaiserslautern
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
 * Authors:
 *  - SciComp, TU Kaiserslautern:
 *     Max Sagebaum
 *     Tim Albring
 *     Johannes Blühdorn
 */

#pragma once

#include "expressionInterface.hpp"

/**
 * @brief Global namespace for CoDiPack - Code Differentiation Package
 */
namespace codi {

  template<typename A> struct ExpressionTraits;

  /**
   * @brief Expression implementation for binary operation with two active variables.
   *
   * @tparam Real The real type used in the active types.
   * @tparam A The expression for the first argument of the function.
   * @tparam B The expression for the second argument of the function.
   * @tparam Impl Implementation of BinaryOpInterface.
   */
  template<typename Real, class A, class B, template<typename> class Impl>
  struct BinaryOp11: public Expression<Real, BinaryOp11<Real, A, B, Impl> > {
    private:

      /** @brief The first argument of the function. */
      typename TypeTraits<A>::StoreType a_;

      /** @brief The second argument of the function. */
      typename TypeTraits<B>::StoreType b_;

    public:
      /**
       * @brief The passive type used in the origin.
       *
       * If Real is not an ActiveReal this value corresponds to Real,
       * otherwise the PassiveValue from Real is used.
       */
      typedef typename TypeTraits<Real>::PassiveReal PassiveReal;

      /**
       * @brief Stores both arguments of the expression.
       *
       * @param[in] a  First argument of the expression.
       * @param[in] b  Second argument of the expression.
       */
      explicit BinaryOp11(const Expression<Real, A>& a, const Expression<Real, B>& b) :
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
          a_.calcGradient(data, Impl<Real>::gradientA(a_.getValue(), b_.getValue(), getValue()));
          b_.calcGradient(data, Impl<Real>::gradientB(a_.getValue(), b_.getValue(), getValue()));
  #else
          Impl<Real>::derv11(data, a_, b_, getValue());
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
          a_.calcGradient(data, Impl<Real>::gradientA(a_.getValue(), b_.getValue(), getValue()) * multiplier);
          b_.calcGradient(data, Impl<Real>::gradientB(a_.getValue(), b_.getValue(), getValue()) * multiplier);
  #else
          Impl<Real>::derv11M(data, a_, b_, getValue(), multiplier);
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
        return Impl<Real>::primal(a_.getValue(), b_.getValue());
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

        return Impl<Real>::primal(aPrimal, bPrimal);
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
        const Real resPrimal = Impl<Real>::primal(aPrimal, bPrimal);

        const PRIMAL_SEED_TYPE aJac = Impl<Real>::gradientA(aPrimal, bPrimal, resPrimal) * seed;
        const PRIMAL_SEED_TYPE bJac = Impl<Real>::gradientB(aPrimal, bPrimal, resPrimal) * seed;
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
        const Real resPrimal = Impl<Real>::primal(aPrimal, bPrimal);

        const Real aJac = Impl<Real>::gradientA(aPrimal, bPrimal, resPrimal) * seed;
        const Real bJac = Impl<Real>::gradientB(aPrimal, bPrimal, resPrimal) * seed;
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
   * @brief Expression implementation for binary operation, only first variable active.
   *
   * @tparam Real The real type used in the active types.
   * @tparam A The expression for the first argument of the function
   * @tparam Impl Implementation of BinaryOpInterface.
   */
  template<typename Real, class A, template<typename> class Impl>
  struct BinaryOp10: public Expression<Real, BinaryOp10<Real, A, Impl> > {
    private:
      /** @brief The type for the passive values. */
      typedef typename TypeTraits<Real>::PassiveReal PassiveReal;

      /** @brief Active first argument of the function */
      typename TypeTraits<A>::StoreType a_;

      /** @brief Passive second argument of the function */
      const PassiveReal b_;

    public:

      /**
       * @brief Stores both arguments of the expression.
       *
       * @param[in] a  First argument of the expression.
       * @param[in] b  Second argument of the expression.
       */
      explicit BinaryOp10(const Expression<Real, A>& a, const PassiveReal& b) :
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
          a_.calcGradient(data, Impl<Real>::gradientA(a_.getValue(), b_, getValue()));
  #else
          Impl<Real>::derv10(data, a_, b_, getValue());
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
          a_.calcGradient(data, Impl<Real>::gradientA(a_.getValue(), b_, getValue()) * multiplier);
  #else
          Impl<Real>::derv10M(data, a_, b_, getValue(), multiplier);
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
        return Impl<Real>::primal(a_.getValue(), b_);
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

        return Impl<Real>::primal(aPrimal, bPrimal);
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
        const Real resPrimal = Impl<Real>::primal(aPrimal, bPrimal);

        const PRIMAL_SEED_TYPE aJac = Impl<Real>::gradientA(aPrimal, bPrimal, resPrimal) * seed;
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
        const Real resPrimal = Impl<Real>::primal(aPrimal, bPrimal);

        const Real aJac = Impl<Real>::gradientA(aPrimal, bPrimal, resPrimal) * seed;
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
   * @brief Expression implementation for binary operation, only second variable active.
   *
   * @tparam Real The real type used in the active types.
   * @tparam B The expression for the second argument of the function.
   * @tparam Impl Implementation of BinaryOpInterface.
   */
  template<typename Real, class B, template<typename> class Impl>
  struct BinaryOp01 : public Expression<Real, BinaryOp01<Real, B, Impl> > {
    private:

      /** @brief The type for the passive values. */
      typedef typename TypeTraits<Real>::PassiveReal PassiveReal;

      /** @brief Passive first argument of the function */
      const PassiveReal a_;

      /** @brief Active second argument of the function */
      typename TypeTraits<B>::StoreType b_;
    public:

      /**
       * @brief Stores both arguments of the expression.
       *
       * @param[in] a  First argument of the expression.
       * @param[in] b  Second argument of the expression.
       */
      explicit BinaryOp01(const PassiveReal& a, const Expression<Real, B>& b) :
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
          b_.calcGradient(data, Impl<Real>::gradientB(a_, b_.getValue(), getValue()));
  #else
          Impl<Real>::derv01(data, a_, b_, getValue());
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
          b_.calcGradient(data, Impl<Real>::gradientB(a_, b_.getValue(), getValue()) * multiplier);
  #else
          Impl<Real>::derv01M(data, a_, b_, getValue(), multiplier);
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
        return Impl<Real>::primal(a_, b_.getValue());
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

        return Impl<Real>::primal(aPrimal, bPrimal);
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
        const Real resPrimal = Impl<Real>::primal(aPrimal, bPrimal);

        const PRIMAL_SEED_TYPE bJac = Impl<Real>::gradientB(aPrimal, bPrimal, resPrimal) * seed;
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
        const Real resPrimal = Impl<Real>::primal(aPrimal, bPrimal);

        const Real bJac = Impl<Real>::gradientB(aPrimal, bPrimal, resPrimal) * seed;
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
   * @brief Specialization of the TypeTraits for the binary operation type.
   *
   * @tparam Real  The floating point value of the active real.
   * @tparam    A  The type of the first argument of the binary operation.
   * @tparam    B  The type of the second argument of the binary operation.
   * @tparam Impl  Implementation of BinaryOpInterface.
   */
  template<typename RealType, typename A, typename B, template<typename> class Impl>
  class TypeTraits< BinaryOp11<RealType, A, B, Impl> > {
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
       * @brief Expressions are temporary and therefore stored by value.
       */
      typedef const BinaryOp11<RealType, A, B, Impl> StoreType;

      /** @brief The maximum derivative order that the active type contains. */
      static const size_t MaxDerivativeOrder = 1 + TypeTraits<Real>::MaxDerivativeOrder;

      /**
       * @brief Get the primal value of the origin of this type.
       * @param[in] t The value from which the primal is extracted.
       * @return The primal value of the origin of this type..
       */
      static const typename TypeTraits<RealType>::PassiveReal getBaseValue(const BinaryOp11<RealType, A, B, Impl>& t) {
        return TypeTraits<RealType>::getBaseValue(t.getValue());
      }
  };

  /**
   * @brief Specialization of the TypeTraits for the binary operation type.
   *
   * @tparam Real  The floating point value of the active real.
   * @tparam    A  The type of the first argument of the binary operation.
   * @tparam Impl  Implementation of BinaryOpInterface.
   */
  template<typename RealType, typename A, template<typename> class Impl>
  class TypeTraits< BinaryOp10<RealType, A, Impl> > {
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
       * @brief Expressions are temporary and therefore stored by value.
       */
      typedef const BinaryOp10<RealType, A, Impl> StoreType;

      /** @brief The maximum derivative order that the active type contains. */
      static const size_t MaxDerivativeOrder = 1 + TypeTraits<Real>::MaxDerivativeOrder;

      /**
       * @brief Get the primal value of the origin of this type.
       * @param[in] t The value from which the primal is extracted.
       * @return The primal value of the origin of this type..
       */
      static const typename TypeTraits<RealType>::PassiveReal getBaseValue(const BinaryOp10<RealType, A, Impl>& t) {
        return TypeTraits<RealType>::getBaseValue(t.getValue());
      }
  };

  /**
   * @brief Specialization of the TypeTraits for the binary operation type.
   *
   * @tparam Real  The floating point value of the active real.
   * @tparam    B  The type of the second argument of the binary operation.
   * @tparam Impl  Implementation of BinaryOpInterface.
   */
  template<typename RealType, typename B, template<typename> class Impl>
  class TypeTraits< BinaryOp01<RealType, B, Impl> > {
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
       * @brief Expressions are temporary and therefore stored by value.
       */
      typedef const BinaryOp01<RealType, B, Impl> StoreType;

      /** @brief The maximum derivative order that the active type contains. */
      static const size_t MaxDerivativeOrder = 1 + TypeTraits<Real>::MaxDerivativeOrder;

      /**
       * @brief Get the primal value of the origin of this type.
       * @param[in] t The value from which the primal is extracted.
       * @return The primal value of the origin of this type..
       */
      static const typename TypeTraits<RealType>::PassiveReal getBaseValue(const BinaryOp01<RealType, B, Impl>& t) {
        return TypeTraits<RealType>::getBaseValue(t.getValue());
      }
  };
}
