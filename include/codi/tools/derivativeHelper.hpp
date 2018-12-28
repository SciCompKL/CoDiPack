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

#include <type_traits>

#include "../configure.h"
#include "binomial.hpp"
#include "../typeTraits.hpp"

/**
 * @brief Global namespace for CoDiPack - Code Differentiation Package
 */
namespace codi {

  namespace DerivativeHelperTemplates {

    /**
     * @brief Creates a reference type from the provided type, that correctly defines the constant modified.
     *
     * The type is defined in type.
     *
     * @tparam     Type  The type for which the reference is created.
     * @tparam constant  Defines if the reference should be to a constant type or not.
     */
    template<typename Type, bool constant>
    struct RefType {

      /** @brief The reference to the specified type */
      typedef typename std::conditional<constant, const Type&, Type&>::type type;
    };

    /**
     * @brief Helper class for the selection of the return type.
     *
     * @tparam  Real  An ActiveReal that has uses the same type for the primal and derivative values.
     * @tparam depth  The selection depth for the return type.
     */
    template<typename Real, size_t depth>
    struct SelectReturnType {

      /** @brief Select the return type in a recursion that walks the primal steps in the value graph. */
      typedef typename SelectReturnType<typename Real::Real, depth -1>::ReturnType ReturnType;
    };

    /**
     * @brief Specialization for the termination of a return type.
     *
     * The termination defines the current node type as the return type.
     *
     * @tparam  Real  An ActiveReal that has uses the same type for the primal and derivative values.
     */
    template<typename Real>
    struct SelectReturnType<Real, 0> {

      /** @brief The current node type is the return type. */
      typedef Real ReturnType;
    };

    /**
     * @brief The walker for the graph of the AD values.
     *
     * The class will recursively walk through the value graph of the AD value.
     *
     * @tparam      Real  The AD type for which the selection should be done.
     *                    Needs to be an ActiveReal type that uses for the primal ad derivative values the same type.
     * @tparam  constant  Defines if the references should be to constant types.
     * @tparam recursion  The number of steps in the graph.
     */
    template<typename Real, bool constant, size_t recursion>
    struct DerivativeSelector {

      /**  @brief The return type for the selection process.  */
      typedef typename SelectReturnType<Real, recursion>::ReturnType ReturnType;

      /** @brief Definition of the ReturnType type as a reference. */
      typedef typename RefType<ReturnType, constant>::type ReturnTypeRef;

      /** @brief Definition of the Real type as an argument. */
      typedef typename RefType<Real, constant>::type RealRef;

      /**
       * @brief Get the specified derivative from the value.
       *
       * The method walks recursively through the graph of the AD type and return the
       * value that is selected with order and l.
       *
       * @param[in] value  The value from which the derivatives are selected.
       * @param[in] order  The order of the derivative that should be selected.
       * @param[in]     l  The index of the derivative that should be returned from the specified order.
       *
       * @return The derivative value that is specified through order and l.
       */
      static ReturnTypeRef select(RealRef value, size_t order, size_t l) {
        size_t lowerDerivatives = binomial(recursion - 1, order);

        if(lowerDerivatives <= l) {
          // take upper branch
          return DerivativeSelector<typename Real::GradientValue, constant, recursion - 1>::select(value.gradient(), order - 1, l - lowerDerivatives);
        } else {
          // take lower branch
          return DerivativeSelector<typename Real::Real, constant, recursion - 1>::select(value.value(), order, l);
        }
      }
    };

    /**
     * @brief The specialization for the termination of the recursion for the walker.
     *
     * This specialization terminates the walk through the value graph.
     *
     * @tparam      Real  The AD type for which the selection should be done.
     *                    Needs to be an ActiveReal type that uses for the primal ad derivative values the same type.
     * @tparam  constant  Defines if the references should be to constant types.
     */
    template<typename Real, bool constant>
    struct DerivativeSelector<Real, constant, 0> {

      /** @brief Definition of the Real type as an argument. */
      typedef typename RefType<Real, constant>::type RealRef;

      /**
       * @brief The termination just return the current selected node.
       *
       * The method walks recursively through the graph of the AD type and return the
       * value that is selected with order and l.
       *
       * @param[in] value  The value from which the derivatives are selected.
       * @param[in] order  The order of the derivative that should be selected.
       * @param[in]     l  The index of the derivative that should be returned from the specified order.
       *
       * @return The derivative value that is specified through order and l.
       */
      static RealRef select(RealRef value, size_t order, size_t l) {
        codiAssert(order == 0); // The order is reduced every time a derivative path is chosen so the order must be 0
        codiAssert(l == 0); // The selection index is adapted every time a path is chosen so the index must be 0

        return value;
      }
    };

    /**
     * @brief The branch selector for the walk through the value graph of an AD type.
     *
     * The branch is selected via lower template argument and specializations for both cases are provided.
     * The lower == true specialization selects the primal value part of the AD type.
     * The lower == true (upper) specialization selects the gradient value part of the AD type.
     *
     * The branch class provides also other values that are required for the walk through the value graph of the AD type.
     *
     * @tparam             Real  The AD type for which the selection should be done. Needs to be an ActiveReal type.
     * @tparam         constant  Defines if the references should be to constant types.
     * @tparam            lower  The selector for the lower or upper branch of the graph. The class is specialized for
     *                           this parameter.
     * @tparam            order  The derivative order that the user wants to select.
     * @tparam                l  The number of the derivative that the user wants to select.
     * @tparam lowerDerivatives  The number of the derivatives, that are available in the lower branch.
     */
    template<typename Real, bool constant, bool lower, size_t order, size_t l, size_t lowerDerivatives>
    struct BranchSelector;

    /**
     * @brief The specialization of the branch selector for lower part of the graph.
     *
     * The specialization selects the lower == true part of the graph.
     *
     * @tparam             Real  The AD type for which the selection should be done. Needs to be an ActiveReal type.
     * @tparam         constant  Defines if the references should be to constant types.
     * @tparam            order  The derivative order that the user wants to select.
     * @tparam                l  The number of the derivative that the user wants to select.
     * @tparam lowerDerivatives  The number of the derivatives, that are available in the lower branch.
     */
    template<typename Real, bool constant, size_t order, size_t l, size_t lowerDerivatives>
    struct BranchSelector<Real, constant, true, order, l, lowerDerivatives> {
      /**
       * @brief The derivative order of the selected branch.
       *
       * It stays the same because we walk along the primal value path.
       */
      static constexpr size_t branchOrder = order;

      /**
       * @brief The number of the derivative that the user wants to select for the branch.
       *
       * It stays the same because we walk along the primal value path.
       */
      static constexpr size_t branchL = l;

      /**
       * @brief The type of the selected branch.
       *
       * Because we walk along the primal value path it is the type of the primal value.
       */
      typedef typename Real::Real Type;

      /** @brief Reference type for the return value. */
      typedef typename RefType<Type, constant>::type TypeRef;
      /** @brief Reference type for the argument value. */
      typedef typename RefType<Real, constant>::type RealRef;

      /**
       * @brief Select the correct branch from the given value.
       *
       * The method will select the primal value of the provided argument.
       *
       * @param[in] value  The value from which the branch is selected.
       */
      static TypeRef select(RealRef value) {
        return value.value();
      }
    };

    /**
     * @brief The specialization of the branch selector for upper part of the graph.
     *
     * The specialization selects the upper (lower == false) part of the graph.
     *
     * @tparam             Real  The AD type for which the selection should be done. Needs to be an ActiveReal type.
     * @tparam         constant  Defines if the references should be to constant types.
     * @tparam            order  The derivative order that the user wants to select.
     * @tparam                l  The number of the derivative that the user wants to select.
     * @tparam lowerDerivatives  The number of the derivatives, that are available in the lower branch.
     */
    template<typename Real, bool constant, size_t order, size_t l, size_t lowerDerivatives>
    struct BranchSelector<Real, constant, false, order, l, lowerDerivatives> {
      /**
       * @brief The derivative order of the selected branch.
       *
       * The order is decrement because we walk along the derivative path and therefore every value in the path has a higher derivative.
       */
      static constexpr size_t branchOrder = order - 1;

      /**
       * @brief The number of the derivative that the user wants to select for the branch.
       *
       * The value is decremented by the number of possible derivatives in the lower branch because these derivatives are ignored by the selection of the upper branch.
       */
      static constexpr size_t branchL = l - lowerDerivatives;

      /**
       * @brief The type of the selected branch.
       *
       * Because we walk along the gradient value path it is the type of the gradient value.
       */
      typedef typename Real::GradientValue Type;

      /** @brief Reference type for the return value. */
      typedef typename std::conditional<constant, const Type&, Type&>::type TypeRef;
      /** @brief Reference type for the argument value. */
      typedef typename std::conditional<constant, const Real&, Real&>::type RealRef;

      /**
       * @brief Select the correct branch from the given value.
       *
       * The method will select the gradient value of the provided argument.
       *
       * @param[in] value  The value from which the branch is selected.
       */
      static TypeRef select(RealRef value) {
        return value.gradient();
      }
    };

    /**
     * @brief The walker for the graph of the AD values
     *
     * The class will recursively walk through the value graph of the AD value.
     *
     * @tparam      Real  The AD type for which the selection should be done.
     *                    Needs to be an ActiveReal type for higher order derivatives.
     *                    For the last step (recursion == 0) it can be any type.
     * @tparam  constant  Defines if the references should be to constant types.
     * @tparam recursion  The number of steps in the graph.
     * @tparam     order  The derivative order that the user wants to select.
     * @tparam         l  The number of the derivative that the user wants to select.
     */
    template<typename Real, bool constant, size_t recursion, size_t order, size_t l>
    struct DerivativeSelectorTemplate {

      /** @brief Compute the number of possible derivatives in the lower branch. */
      static constexpr size_t lowerDerivatives = Binomial<recursion - 1, order>::value;
      /** @brief Check if the lower branch should be selected */
      static constexpr   bool lowerBranch = lowerDerivatives > l;

      /** @brief Choose the selected branch. */
      typedef BranchSelector<Real, constant, lowerBranch, order, l, lowerDerivatives> Branch;

      /** @brief define the next step in the graph walk */
      typedef DerivativeSelectorTemplate<typename Branch::Type, constant, recursion - 1, Branch::branchOrder, Branch::branchL> NextStep;
      /** @brief Recursively define the return type of the selected value in the graph */
      typedef typename NextStep::ReturnTypeRef ReturnTypeRef;


      /** @brief Reference type for the argument value. */
      typedef typename RefType<Real, constant>::type RealRef;

      /**
       * @brief Select the correct branch from the given value.
       *
       * @param[in] value The current value for which the graph is walked.
       *
       * @return The selected value in the graph.
       */
      static ReturnTypeRef select(RealRef value) {
        return NextStep::select(Branch::select(value));
      }
    };

    /**
     * @brief Specialization for the final step in the graph walk.
     *
     * In the final step the current value is returned.
     *
     * @tparam     Real  The AD type for which the selection should be done.
     *                   Needs to be an ActiveReal type for higher order derivatives.
     *                   For the last step (recursion == 0) it can be any type.
     * @tparam constant  Defines if the references should be to constant types.
     * @tparam    order  The derivative order that the user wants to select.
     * @tparam        l  The number of the derivative that the user wants to select.
     */
    template<typename Real, bool constant, size_t order, size_t l>
    struct DerivativeSelectorTemplate<Real, constant, 0, order, l> {

      /** @brief The final return type is the current type of the graph node. */
      typedef typename std::conditional<constant, const Real&, Real&>::type ReturnTypeRef;

      /**
       * @brief The final step just return the current value of the node.
       *
       * @param[in] value The current value for which the graph is walked.
       *
       * @return The selected value in the graph.
       */
      static ReturnTypeRef select(ReturnTypeRef value) {
        return value;
      }
    };

    /**
     * @brief Compile time loop for the SetDerivatives method
     *
     * The loop will decrement the parameter l until it reaches zero.
     *
     * It sets the derivative on for l-1, l-2, ..., 0.
     *
     * @tparam      Real  The AD type for which the selection should be done.
     *                    Needs to be an ActiveReal type.
     * @tparam      Type  The type of the value that is set.
     * @tparam recursion  The number of steps in the graph.
     * @tparam     order  The derivative order that the user wants to select.
     * @tparam         l  The number of the derivative that the user wants to select.
     */
    template<typename Real, typename Type, size_t recursion, size_t order, size_t l>
    struct SetDerivativesLoop {

      /**
       * @brief perform the loop in a recursive fashion
       *
       * @param[in,out]      value  The value on which the derivatives are set.
       * @param[in]     derivative  The derivative that is set in the value graph.
       */
      static void doLoop(Real& value, const Type& derivative) {
        DerivativeSelectorTemplate<Real, false, recursion, order, l - 1>::select(value) = derivative;

        SetDerivativesLoop<Real, Type, recursion, order, l - 1>::doLoop(value, derivative);
      }
    };

    /**
     * @brief Specialization for the termination of the compile time loop for the SetDerivatives method
     *
     * @tparam      Real  The AD type for which the selection should be done.
     *                    Needs to be an ActiveReal type.
     * @tparam      Type  The type of the value that is set.
     * @tparam recursion  The number of steps in the graph.
     * @tparam     order  The derivative order that the user wants to select.
     */
    template<typename Real, typename Type, size_t recursion, size_t order>
    struct SetDerivativesLoop<Real, Type, recursion, order, 0> {

      /**
       * @brief Terminates the loop.
       *
       * @param[in,out]      value  The value on which the derivatives are set.
       * @param[in]     derivative  The derivative that is set in the value graph.
       */
      static void doLoop(Real& value, const Type& derivative) {
        CODI_UNUSED(value);
        CODI_UNUSED(derivative);

        // termination of the loop
      }
    };
  }

  /**
   * @brief A helper class for the convenient selection of gradient data of higher order AD types.
   *
   * The algorithm behind the selection walks along the value path of the derivatives. A higher order
   * derivative, that is combined via the CoDiPack types, has \f$2^n\f$ possible derivative values
   * (including the primal value) with \f$n\f$ the number of nested types. If a second and third order type is
   * constructed via
   * \code{.cpp}
   *  typedef RealForwardGen<RealForward> t2s;
   *  typedef RealForwardGen<t2s>         t3s;
   * \endcode
   * then the second order type t2s has 4 possible values \f$(n = 2)\f$ and the third order type has 8
   * possible values \f$(n = 3)\f$. The number of derivatives per derivative order can be computed via
   * the binomial coefficient \f$(n over k)\f$ where \f$k\f$ is the derivative order. For a second order type
   * this yields 1 derivative of zero order (the primal value) two derivatives of the first order and one
   * derivative of the second order. For a third order type this is 0:1, 1:3, 2:3, 3:1.
   *
   * The algorithm in the class will now select the appropriate derivative. The user provides the values
   * k and the number of the derivative he wants to select i.e. \f$0, ... (n over k) - 1\f$.
   *
   * The algorithm walks the graph of the CoDiPack types and selects if the primal value or the derivative
   * value should be selected. For the third order type the graph looks like:
   *
   *  t3s  t2s  t1s  double | order  index
   *                        |
   *               ,---o    |  3     0
   *              /         |
   *            ,o-----o    |  2     2
   *           /            |
   *          /    ,---o    |  2     1
   *         /    /         |
   *        o----o-----o    |  1     2
   *       /                |
   *      /        ,---o    |  2     0
   *     /        /         |
   *    /       ,o-----o    |  1     1
   *   /       /            |
   *   |      /    ,---o    |  1     0
   *   |     /    /         |
   *   o----o----o-----o    |  0     0
   *
   * The lower branch always stand for the primal value and the upper branch always stands for the derivative
   * value.
   * The two columns at the end show the derivative order under the column 'order' and the index in that order
   * class under the column 'index'. It can be seen that the different derivative values of the same order are
   * not continuously ordered in the graph.
   *
   * The class provides  methods for the runtime selection of the derivatives. If theses methods are
   * used, then the AD types need to be defined such that all primal and derivative values have the same type.
   * If this is not the case, then the compiler will show errors, that it can not convert a value.
   *
   * The class provides methods for the compile time selection of the derivatives, too. These methods do not have
   * the restriction, that all the primal and derivative types need to have the same type. On the other hand
   * all compile time restrictions apply to the parameters of the templates. That is, they need to be compile
   * time constants.
   * Also the setDerivatives method which sets all derivatives of one order may not be used if different primal
   * and gradient types are used. The provided objects needs to be convertible into all possible types, that
   * are used at the termination points of the recursion.
   *
   * @tparam Real  The AD type for which the derivatives are selected. The type needs to be an ActiveReal.
   */
  template<typename Real>
  struct DerivativeHelper {

    /**
     * @brief Selection of a derivative component from the given value.
     *
     * This methods requires the real type to use the same value for the primal values and the
     * derivative values.
     *
     * @param[in] value  The value from which the derivative is selected.
     * @param[in] order  The order of the derivative that should be selected.
     * @param[in]     l  The number of the derivative that should be selected in the given order.
     *
     * @return The reference to the selected derivative value.
     *
     * @tparam depth  The selection depth of the algorithm. The default value is the maximum derivative order.
     *                If a smaller depth is provided, subgraphs can be chosen.
     */
    template<size_t depth = TypeTraits<Real>::MaxDerivativeOrder>
    static typename DerivativeHelperTemplates::DerivativeSelector<Real, false, depth>::ReturnTypeRef derivative(Real& value, size_t order, size_t l) {

      if(order > TypeTraits<Real>::MaxDerivativeOrder) {
        CODI_EXCEPTION("The derivative order must be smaller or equal than the maximum provided derivative. order: %d, max derivative: %d.", order, TypeTraits<Real>::MaxDerivativeOrder);
      }

      size_t numberDerivatives = binomial(TypeTraits<Real>::MaxDerivativeOrder, order);
      if(l >= numberDerivatives) {
        CODI_EXCEPTION("The selected derivative must be smaller than the maximum number of derivatives. selected: %d, number derivatives: %d.", l, numberDerivatives);
      }

      return DerivativeHelperTemplates::DerivativeSelector<Real, false, depth>::select(value, order , l);
    }

    /**
     * @brief Selection of a derivative component from the given value.
     *
     * This methods requires the real type to use the same value for the primal values and the
     * derivative values.
     *
     * @param[in] value  The value from which the derivative is selected.
     * @param[in] order  The order of the derivative that should be selected.
     * @param[in]     l  The number of the derivative that should be selected in the given order.
     *
     * @return The reference to the selected derivative value.
     *
     * @tparam depth  The selection depth of the algorithm. The default value is the maximum derivative order.
     *                If a smaller depth is provided, subgraphs can be chosen.
     */
    template<size_t depth = TypeTraits<Real>::MaxDerivativeOrder>
    static typename DerivativeHelperTemplates::DerivativeSelector<Real, true, depth>::ReturnTypeRef derivative(const Real& value, size_t order, size_t l) {

      if(order > TypeTraits<Real>::MaxDerivativeOrder) {
        CODI_EXCEPTION("The derivative order must be smaller or equal than the maximum provided derivative. order: %d, max derivative: %d.", order, TypeTraits<Real>::MaxDerivativeOrder);
      }

      size_t numberDerivatives = binomial(TypeTraits<Real>::MaxDerivativeOrder, order);
      if(l >= numberDerivatives) {
        CODI_EXCEPTION("The selected derivative must be smaller than the maximum number of derivatives. selected: %d, number derivatives: %d.", l, numberDerivatives);
      }

      return DerivativeHelperTemplates::DerivativeSelector<Real, true, depth>::select(value, order , l);
    }

    /**
     * @brief Set the derivative for all values of a given order.
     *
     * The method iterates over the value from 0 to (depth over order)-1.
     *
     * This method may cause compilation errors if mixed primal and derivative types are used
     * in the ActiveReal type. See the DerivativeHelper documentation for details.
     *
     * @param[in,out]      value  The value on which the derivatives are set.
     * @param[in]          order  The order of the derivative that should be selected.
     * @param[in]     derivative  The value that is set on the derivatives.
     *
     * @tparam depth  The selection depth of the algorithm. The default value is the maximum derivative order.
     *                If a smaller depth is provided, subgraphs can be chosen.
     */
    template<size_t depth = TypeTraits<Real>::MaxDerivativeOrder>
    static void setDerivatives(Real& value, size_t order, const typename DerivativeHelperTemplates::DerivativeSelector<Real, false, depth>::ReturnType& derivative) {

      if(order > TypeTraits<Real>::MaxDerivativeOrder) {
        CODI_EXCEPTION("The derivative order must be smaller or equal than the maximum provided derivative. order: %d, max derivative: %d.", order, TypeTraits<Real>::MaxDerivativeOrder);
      }

      size_t numberDerivatives = binomial(TypeTraits<Real>::MaxDerivativeOrder, order);
      for(size_t i = 0; i < numberDerivatives; ++i) {
        DerivativeHelperTemplates::DerivativeSelector<Real, false, depth>::select(value, order , i) = derivative;
      }
    }

    /**
     * @brief Set the derivative for all values of a given order for the tape recording.
     *
     * The method is the same as setDerivatives but only for the directions that are defined
     * during the tape recording. The method assumes that the provided type as based on a reverse
     * tape implementation and that all nested types are forward tape implementations.
     *
     * This method may cause compilation errors if mixed primal and derivative types are used
     * in the ActiveReal type. See the DerivativeHelper documentation for details.
     *
     * @param[in,out]      value  The value on which the derivatives are set.
     * @param[in]          order  The order of the derivative that should be selected.
     * @param[in]     derivative  The value that is set on the derivatives.
     *
     * @tparam depth  The selection depth of the algorithm. The default value is the maximum derivative order.
     *                If a smaller depth is provided, subgraphs can be chosen.
     */
    template<size_t depth = TypeTraits<Real>::MaxDerivativeOrder>
    static void setDerivativesForward(Real& value, size_t order, const typename DerivativeHelperTemplates::DerivativeSelector<Real, false, depth>::ReturnType& derivative) {
      if(order > TypeTraits<Real>::MaxDerivativeOrder - 1) {
        CODI_EXCEPTION("The derivative order must be smaller or equal than the maximum provided forward derivative. order: %d, max forward derivative: %d.", order, TypeTraits<Real>::MaxDerivativeOrder - 1);
      }

      DerivativeHelper<typename Real::Real>::template setDerivatives<depth - 1>(value.value(), order, derivative);
    }

    /**
     * @brief Set the derivative for all values of a given order for the tape reverse interpretation.
     *
     * The method is the same as setDerivatives but only for the directions that are not yet defined
     * during the tape recording. It sets all missing directions for the reverse interpretation.
     * The method assumes that the provided type as based on a reverse tape implementation and that
     * all nested types are forward tape implementations.
     *
     * This method may cause compilation errors if mixed primal and derivative types are used
     * in the ActiveReal type. See the DerivativeHelper documentation for details.
     *
     * @param[in,out]      value  The value on which the derivatives are set.
     * @param[in]          order  The order of the derivative that should be selected.
     * @param[in]     derivative  The value that is set on the derivatives.
     *
     * @tparam depth  The selection depth of the algorithm. The default value is the maximum derivative order.
     *                If a smaller depth is provided, subgraphs can be chosen.
     */
    template<size_t depth = TypeTraits<Real>::MaxDerivativeOrder>
    static void setDerivativesReverse(Real& value, size_t order, const typename DerivativeHelperTemplates::DerivativeSelector<Real, false, depth>::ReturnType& derivative) {

      if(order > TypeTraits<Real>::MaxDerivativeOrder) {
        CODI_EXCEPTION("The derivative order must be smaller or equal than the maximum provided reverse derivative. order: %d, max reverse derivative: %d.", order, TypeTraits<Real>::MaxDerivativeOrder);
      }
      if(0 == order) {
        CODI_EXCEPTION("The derivative order must be at least one for reverse derivatives. order: %d.", order);
      }

      DerivativeHelper<typename Real::GradientValue>::template setDerivatives<depth - 1>(value.gradient(), order - 1, derivative);
    }

    /**
     * @brief Selection of a derivative component from the given value.
     *
     * @param[in] value  The value from which the derivative is selected.
     *
     * @return The reference to the selected derivative value.
     *
     * @tparam order  The order of the derivative that should be selected.
     * @tparam     l  The number of the derivative that should be selected in the given order.
     * @tparam depth  The selection depth of the algorithm. The default value is the maximum derivative order.
     *                If a smaller depth is provided, subgraphs can be chosen.
     */
    template<size_t order, size_t l, size_t depth = TypeTraits<Real>::MaxDerivativeOrder>
    static typename DerivativeHelperTemplates::DerivativeSelectorTemplate<Real, false, depth, order, l>::ReturnTypeRef derivative(Real& value) {
      return DerivativeHelperTemplates::DerivativeSelectorTemplate<Real, false, depth, order, l>::select(value);
    }

    /**
     * @brief Selection of a derivative component from the given value.
     *
     * @param[in] value  The value from which the derivative is selected.
     *
     * @return The constant reference to the selected derivative value.
     *
     * @tparam order  The order of the derivative that should be selected.
     * @tparam     l  The number of the derivative that should be selected in the given order.
     * @tparam depth  The selection depth of the algorithm. The default value is the maximum derivative order.
     *                If a smaller depth is provided, subgraphs can be chosen.
     */
    template<size_t order, size_t l, size_t depth = TypeTraits<Real>::MaxDerivativeOrder>
    static typename DerivativeHelperTemplates::DerivativeSelectorTemplate<Real, true, depth, order, l>::ReturnTypeRef derivative(const Real& value) {
      return DerivativeHelperTemplates::DerivativeSelectorTemplate<Real, true, depth, order, l>::select(value);
    }

    /**
     * @brief Set the derivative for all values of a given order.
     *
     * The method iterates over the value from (depth over order)-1 to 0.
     *
     * This method may cause compilation errors if mixed primal and derivative types are used
     * in the ActiveReal type. See the DerivativeHelper documentation for details.
     *
     * @param[in,out]      value  The value on which the derivatives are set.
     * @param[in]     derivative  The value that is set on the derivatives.
     *
     * @tparam order  The order of the derivative that should be selected.
     * @tparam  Type  The type of the value that is set to the derivatives.
     * @tparam depth  The selection depth of the algorithm. The default value is the maximum derivative order.
     *                If a smaller depth is provided, subgraphs can be chosen.
     */
    template<size_t order, typename Type, size_t depth = TypeTraits<Real>::MaxDerivativeOrder>
    static void setDerivatives(Real& value, const Type& derivative) {
      DerivativeHelperTemplates::SetDerivativesLoop<Real, Type, depth, order, Binomial<depth, order>::value >::doLoop(value, derivative);
    }

    /**
     * @brief Set the derivative for all values of a given order for the tape recording.
     *
     * The method is the same as setDerivatives but only for the directions that are defined
     * during the tape recording. The method assumes that the provided type as based on a reverse
     * tape implementation and that all nested types are forward tape implementations.
     *
     * This method may cause compilation errors if mixed primal and derivative types are used
     * in the ActiveReal type. See the DerivativeHelper documentation for details.
     *
     * @param[in,out]      value  The value on which the derivatives are set.
     * @param[in]     derivative  The value that is set on the derivatives.
     *
     * @tparam order  The order of the derivative that should be selected.
     * @tparam  Type  The type of the value that is set to the derivatives.
     * @tparam depth  The selection depth of the algorithm. The default value is the maximum derivative order.
     *                If a smaller depth is provided, subgraphs can be chosen.
     */
    template<size_t order, typename Type, size_t depth = TypeTraits<Real>::MaxDerivativeOrder>
    static void setDerivativesForward(Real& value, const Type& derivative) {
      DerivativeHelper<typename Real::Real>::template setDerivatives<order, Type, depth - 1>(value.value(), derivative);
    }

    /**
     * @brief Set the derivative for all values of a given order for the tape reverse interpretation.
     *
     * The method is the same as setDerivatives but only for the directions that are not yet defined
     * during the tape recording. It sets all missing directions for the reverse interpretation.
     * The method assumes that the provided type as based on a reverse tape implementation and that
     * all nested types are forward tape implementations.
     *
     * This method may cause compilation errors if mixed primal and derivative types are used
     * in the ActiveReal type. See the DerivativeHelper documentation for details.
     *
     * @param[in,out]      value  The value on which the derivatives are set.
     * @param[in]     derivative  The value that is set on the derivatives.
     *
     * @tparam order  The order of the derivative that should be selected.
     * @tparam  Type  The type of the value that is set to the derivatives.
     * @tparam depth  The selection depth of the algorithm. The default value is the maximum derivative order.
     *                If a smaller depth is provided, subgraphs can be chosen.
     */
    template<size_t order, typename Type, size_t depth = TypeTraits<Real>::MaxDerivativeOrder>
    static void setDerivativesReverse(Real& value, const Type& derivative) {
      DerivativeHelper<typename Real::GradientValue>::template setDerivatives<order - 1, Type, depth - 1>(value.gradient(), derivative);
    }
  };
}
