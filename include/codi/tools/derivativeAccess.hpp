/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015-2023 Chair for Scientific Computing (SciComp), University of Kaiserslautern-Landau
 * Homepage: http://www.scicomp.uni-kl.de
 * Contact:  Prof. Nicolas R. Gauger (codi@scicomp.uni-kl.de)
 *
 * Lead developers: Max Sagebaum, Johannes Blühdorn (SciComp, University of Kaiserslautern-Landau)
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
 * For other licensing options please contact us.
 *
 * Authors:
 *  - SciComp, University of Kaiserslautern-Landau:
 *    - Max Sagebaum
 *    - Johannes Blühdorn
 *    - Former members:
 *      - Tim Albring
 */
#pragma once

#include <utility>

#include "../config.h"
#include "../expressions/lhsExpressionInterface.hpp"
#include "../misc/binomial.hpp"
#include "../misc/compileTimeLoop.hpp"
#include "../misc/exceptions.hpp"
#include "../misc/macros.hpp"
#include "../traits/realTraits.hpp"

/** \copydoc codi::Namespace */
namespace codi {

#ifndef DOXYGEN_DISABLE

  /**
   * See DerivativeAccess for details about the number for each derivative order and the iteration order.
   */
  namespace DerivativeAccessImpl {
    CODI_INLINE size_t constexpr maximumDerivatives(size_t selectionDepth, size_t order) {
      return binomial(selectionDepth, order);
    }

    CODI_INLINE size_t constexpr maximumDerivativesPrimalBranch(size_t selectionDepth, size_t order) {
      return binomial(selectionDepth - 1, order);
    }

    CODI_INLINE size_t constexpr isPrimalBranch(size_t selectionDepth, size_t order, size_t l) {
      return l < maximumDerivativesPrimalBranch(selectionDepth, order);
    }

    template<typename Type, size_t selectionDepth, size_t order, size_t l>
    struct CheckCompileTimeValues {
      public:
        CODI_STATIC_ASSERT(selectionDepth <= RealTraits::MaxDerivativeOrder<Type>(),
                           "Selection depth must not be higher than the maximum derivative order.");
        CODI_STATIC_ASSERT(order <= selectionDepth, "Derivative order must not be higher than the selection depth.");
        CODI_STATIC_ASSERT(l < maximumDerivatives(selectionDepth, order),
                           "Selected derivative must not be greater than the number of available derivatives for that"
                           "order.");

        static bool constexpr isValid = true;
    };

    /**
     * @brief Compile time selection of correct derivative.
     *
     * Calls itself recursively.
     *
     * @tparam T_Type CoDiPack type.
     * @tparam constant  If the argument has the constant modifier.
     * @tparam selectionDepth  Maximum derivative order of T_Type.
     * @tparam order  Order of the derivative.
     * @tparam l  Index of the derivative of the chosen order.
     * @tparam primalBranch  Compile time selection of the primal or derivative branch.
     */
    template<typename Type, bool constant, size_t selectionDepth, size_t order, size_t l,
             bool primalBranch = isPrimalBranch(selectionDepth, order, l)>
    struct SelectCompileTime;

    /// \copydoc SelectCompileTime
    template<typename T_Type, bool constant, size_t selectionDepth, size_t order, size_t l>
    struct SelectCompileTime<T_Type, constant, selectionDepth, order, l, true> {
      public:
        using Type = CODI_DD(T_Type, CODI_DEFAULT_LHS_EXPRESSION);

        using Inner = SelectCompileTime<typename Type::Real, constant, selectionDepth - 1, order, l>;
        using ArgType = typename std::conditional<constant, Type const, Type>::type;
        using RType = typename Inner::RType;

        CODI_STATIC_ASSERT(CODI_T(CheckCompileTimeValues<Type, selectionDepth, order, l>::isValid),
                           "Checks inside of type.");

        static RType& select(ArgType& value) {
          return Inner::select(value.value());
        }
    };

    /// \copydoc SelectCompileTime
    template<typename T_Type, bool constant, size_t selectionDepth, size_t order, size_t l>
    struct SelectCompileTime<T_Type, constant, selectionDepth, order, l, false> {
      public:
        using Type = CODI_DD(T_Type, CODI_DEFAULT_LHS_EXPRESSION);

        using Inner = SelectCompileTime<typename Type::Gradient, constant, selectionDepth - 1, order - 1,
                                        l - maximumDerivativesPrimalBranch(selectionDepth, order)>;
        using ArgType = typename std::conditional<constant, Type const, Type>::type;
        using RType = typename Inner::RType;

        CODI_STATIC_ASSERT(CODI_T(CheckCompileTimeValues<Type, selectionDepth, order, l>::isValid),
                           "Checks inside of type.");

        static RType& select(ArgType& value) {
          return Inner::select(value.gradient());
        }
    };

    /// Terminator of selection recursion.
    template<typename Type, bool constant>
    struct SelectCompileTime<Type, constant, 0, 0, 0, true> {
      public:
        using ArgType = typename std::conditional<constant, Type const, Type>::type;
        using RType = ArgType;

        static RType& select(ArgType& value) {
          return value;
        }
    };

    /**
     * @brief Run time selection of a specific derivative.
     *
     * Calls itself recursively.
     *
     * @tparam T_Type CoDiPack type.
     * @tparam constant  If the argument has the constant modifier.
     * @tparam T_selectionDepth  Maximum derivative order of T_Type.
     */
    template<typename T_Type, bool constant, size_t T_selectionDepth>
    struct SelectRunTime {
      public:
        using Type = CODI_DD(T_Type, CODI_DEFAULT_LHS_EXPRESSION);
        static size_t constexpr selectionDepth = CODI_DD(T_selectionDepth, CODI_UNDEFINED_VALUE);

        CODI_STATIC_ASSERT(CODI_T(std::is_same<typename Type::Real, typename Type::Gradient>::value),
                           "CoDiPack type needs to have the same real and gradient value for run time derivative "
                           "selection.");
        CODI_STATIC_ASSERT(selectionDepth <= RealTraits::MaxDerivativeOrder<Type>(),
                           "Selection depth must not be higher than the maximum derivative order");

        using Inner = SelectRunTime<typename Type::Real, constant, selectionDepth - 1>;
        using ArgType = typename std::conditional<constant, Type const, Type>::type;
        using RType = typename Inner::RType;

        static RType& select(ArgType& v, size_t order, size_t l) {
          size_t const maxDerivativesPrimalBranch = binomial(selectionDepth - 1, order);
          if (l < maxDerivativesPrimalBranch) {
            return Inner::select(v.value(), order, l);
          } else {
            return Inner::select(v.gradient(), order - 1, l - maxDerivativesPrimalBranch);
          }
        }
    };

    /// Terminator of selection recursion.
    template<typename T_Type, bool constant>
    struct SelectRunTime<T_Type, constant, 0> {
      public:
        using Type = CODI_DD(T_Type, double);
        using ArgType = typename std::conditional<constant, Type const, Type>::type;
        using RType = ArgType;

        static RType& select(ArgType& v, size_t order, size_t l) {
          CODI_UNUSED(order, l);
          return v;
        }
    };
  }
#endif

  /**
   * @brief A helper class for the access of the various derivatives in higher order AD types.
   *
   * A higher order derivative that is combined via the CoDiPack types has \f$2^n\f$ possible derivative values
   * (including the primal value) where \f$n\f$ is the number of types that are nested. If a second and third order type
   * are constructed via
   * \code{.cpp}
   *  typedef RealForwardGen<RealForward> t2s;
   *  typedef RealForwardGen<t2s>         t3s;
   * \endcode
   * then the second order type t2s has 4 possible values \f$(n = 2)\f$ and the third order type has 8
   * possible values \f$(n = 3)\f$. The number of derivatives per derivative order can be computed via
   * the binomial coefficient \f$\binom{n}{k}\f$ where \f$k\f$ is the derivative order. For a second order type
   * this yields 1 derivative of order zero (the primal value), two derivatives of first order, and one
   * derivative of second order. For a third order type this is 0:1, 1:3, 2:3, 3:1 (order:number).
   *
   * Given the derivative order \f$k\f$ and the index of the derivative within that order
   * \f$l\in\{0, ..., \binom{n}{k} - 1\}\f$, the algorithm in this class selects the corresponding derivative value.
   *
   * The class provides methods for the run time selection of the derivatives. If these methods are
   * used, then the AD types need to be defined such that all primal and derivative values have the same type.
   * If this is not the case, then the compiler will show errors that it cannot convert a value.
   *
   * The class provides methods for the compile time selection of the derivatives, too. These methods do not have
   * the restriction that all the primal and derivative types need to have the same type. On the other hand, the usual
   * compile time restrictions apply to the parameters of the templates, that is, they must be compile time constants.
   * Also, the setDerivatives method which sets all derivatives of one order may not be used if different primal
   * and gradient types are used. The provided objects must be convertible into all possible types that
   * are used at the leaf points of the recursion.
   *
   * The selection algorithm walks along the nodes in the nested classes. From the specified order and the derivative
   * number \f$l\f$ the algorithm checks if the derivative \f$l\f$ is in the lower (value) or upper (gradient) branch of
   * the nested classes. The value is compared against the number of derivatives of that order in the lower branch. If
   * smaller, the lower branch is taken. If larger or equal, the upper branch is taken.
   * For a third order type, the graph looks as follows.
   *
   * \code{.cpp}
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
   * \endcode
   *
   * The two columns at the end show the derivative order under the column 'order' and the index in that order
   * class under the column 'index'. It can be seen that the different derivative values of the same order are
   * not continuously ordered in the graph. For the 2nd order derivatives in the above example, the indices correspond
   * to the following calls.
   *
   * \code{.cpp}
   *   DA::derivative<2,0>(v) == v.value().gradient().gradient();
   *   DA::derivative<2,1>(v) == v.gradient().value().gradient();
   *   DA::derivative<2,2>(v) == v.gradient().gradient().value();
   * \endcode
   *
   * @tparam T_Type  The AD type for which the derivatives are selected. The type must implement LhsExpressionInterface.
   */
  template<typename T_Type>
  struct DerivativeAccess {
    public:

      /// See DerivativeAccess.
      using Type = CODI_DD(T_Type, CODI_DEFAULT_LHS_EXPRESSION);

      /// Helper for the run time selection of derivatives.
      template<bool constant, size_t selectionDepth>
      using SelectRunTime = DerivativeAccessImpl::SelectRunTime<Type, constant, selectionDepth>;

      /// Helper for the compile time selection of derivatives.
      template<bool constant, size_t selectionDepth, size_t order, size_t l>
      using SelectCompileTime = DerivativeAccessImpl::SelectCompileTime<Type, constant, selectionDepth, order, l>;

      /// Run time selection of derivatives. \f$order \in \{0, ..., selectionDepth\}\f$.\f$l \in \{0 ...
      /// (\binom{selectionDepth}{order} - 1)\}\f$.
      template<size_t selectionDepth = RealTraits::MaxDerivativeOrder<Type>()>
      static typename SelectRunTime<true, selectionDepth>::RType const& derivative(Type const& v, size_t order,
                                                                                   size_t l) {
        checkRuntimeSelection<selectionDepth>(order, l);

        return SelectRunTime<true, selectionDepth>::select(v, order, l);
      }

      /// Run time selection of derivatives. \f$order \in \{0, ..., selectionDepth\}\f$. \f$l \in \{0 ...
      /// (\binom{selectionDepth}{order} - 1)\}\f$.
      template<size_t selectionDepth = RealTraits::MaxDerivativeOrder<Type>()>
      static typename SelectRunTime<false, selectionDepth>::RType& derivative(Type& v, size_t order, size_t l) {
        checkRuntimeSelection<selectionDepth>(order, l);

        return SelectRunTime<false, selectionDepth>::select(v, order, l);
      }

      /// Run time set of all derivatives of the same order. \f$order \in \{0, ..., selectionDepth\}\f$.
      template<typename Derivative, size_t selectionDepth = RealTraits::MaxDerivativeOrder<Type>()>
      static void setAllDerivatives(Type& v, size_t order, Derivative const& d) {
        size_t const maxDerivatives = binomial(selectionDepth, order);
        for (size_t i = 0; i < maxDerivatives; i += 1) {
          derivative<selectionDepth>(v, order, i) = d;
        }
      }

      /// Run time set of all derivatives in the primal value part of the same order.
      /// \f$order \in \{0, ..., selectionDepth\}\f$.
      template<typename Derivative, size_t selectionDepth = RealTraits::MaxDerivativeOrder<Type>()>
      static void setAllDerivativesForward(Type& v, size_t order, Derivative const& d) {
        DerivativeAccess<typename Type::Real>::template setAllDerivatives<Derivative, selectionDepth - 1>(v.value(),
                                                                                                          order, d);
      }

      /// Run time set of all derivatives in the gradient part of the same order.
      /// \f$order\in\{0, ..., selectionDepth\}\f$.
      template<typename Derivative, size_t selectionDepth = RealTraits::MaxDerivativeOrder<Type>()>
      static void setAllDerivativesReverse(Type& v, size_t order, Derivative const& d) {
        DerivativeAccess<typename Type::Gradient>::template setAllDerivatives<Derivative, selectionDepth - 1>(
            v.gradient(), order - 1, d);
      }

      /// Compile time selection of derivatives. \f$order \in \{0, ..., selectionDepth\}\f$. \f$l \in \{0 ...
      /// (\binom{selectionDepth}{order} - 1)\}\f$.
      template<size_t order, size_t l, size_t selectionDepth = RealTraits::MaxDerivativeOrder<Type>()>
      static typename SelectCompileTime<true, selectionDepth, order, l>::RType const& derivative(Type const& v) {
        return SelectCompileTime<true, selectionDepth, order, l>::select(v);
      }

      /// Compile time selection of derivatives. \f$order \in \{0, ..., selectionDepth\}\f$. \f$l \in \{0 ...
      /// (\binom{selectionDepth}{order} - 1)\}\f$.
      template<size_t order, size_t l, size_t selectionDepth = RealTraits::MaxDerivativeOrder<Type>()>
      static typename SelectCompileTime<false, selectionDepth, order, l>::RType& derivative(Type& v) {
        return SelectCompileTime<false, selectionDepth, order, l>::select(v);
      }

      /// Compile time set of all derivatives of the same order. \f$order \in \{0, ..., selectionDepth\}\f$.
      template<size_t order, typename Derivative, size_t selectionDepth = RealTraits::MaxDerivativeOrder<Type>()>
      static void setAllDerivatives(Type& v, Derivative const& d) {
        CompileTimeLoop<DerivativeAccessImpl::maximumDerivatives(selectionDepth, order)>::eval(
            CallSetDerivative<order, Derivative, selectionDepth>{}, v, d);
      }

      /// Compile time set of all derivatives in the primal value part of the same order.
      /// \f$order \in \{0, ..., selectionDepth\}\f$.
      template<size_t order, typename Derivative, size_t selectionDepth = RealTraits::MaxDerivativeOrder<Type>()>
      static void setAllDerivativesForward(Type& v, Derivative const& d) {
        DerivativeAccess<typename Type::Real>::template setAllDerivatives<order, Derivative, selectionDepth - 1>(
            v.value(), d);
      }

      /// Compile time set of all derivatives in the gradient part of the same order.
      /// \f$order \in \{0, ..., selectionDepth\}\f$.
      template<size_t order, typename Derivative, size_t selectionDepth = RealTraits::MaxDerivativeOrder<Type>()>
      static void setAllDerivativesReverse(Type& v, Derivative const& d) {
        DerivativeAccess<typename Type::Gradient>::template setAllDerivatives<order - 1, Derivative,
                                                                              selectionDepth - 1>(v.gradient(), d);
      }

    private:

      template<size_t selectionDepth>
      static void checkRuntimeSelection(size_t order, size_t l) {
        if (order > selectionDepth) {
          CODI_EXCEPTION(
              "The derivative order must be smaller or equal than the maximum possible derivative. order: %d, max "
              "derivative: %d.",
              order, selectionDepth);
        }

        size_t numberDerivatives = binomial(selectionDepth, order);
        if (l >= numberDerivatives) {
          CODI_EXCEPTION(
              "The selected derivative must be smaller than the maximum number of derivatives. selected: %d, number "
              "derivatives: %d.",
              l, numberDerivatives);
        }
      }

      template<size_t order, typename Derivative, size_t selectionDepth>
      struct CallSetDerivative {
        public:
          template<size_t pos>
          void operator()(std::integral_constant<size_t, pos>, Type& v, Derivative const& d) {
            DerivativeAccess::derivative<order, pos - 1, selectionDepth>(v) = d;
          }
      };
  };
}
