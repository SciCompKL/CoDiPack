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

#include <type_traits>

#include "../config.h"
#include "../expressions/logic/compileTimeTraversalLogic.hpp"
#include "../misc/macros.hpp"
#include "misc/enableIfHelpers.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  template<typename T_Real, typename T_Impl>
  struct ExpressionInterface;

  template<typename T_Real, typename T_Gradient, typename T_Tape, typename T_Impl>
  struct LhsExpressionInterface;

  template<typename T_Real, template<typename> class T_ConversionOperator>
  struct ConstantExpression;

  template<typename T_Tape>
  struct StaticContextActiveType;

  /// Traits for everything that can be an expression e.g. codi::RealReverse, a + b, etc..
  namespace ExpressionTraits {

    /*******************************************************************************/
    /// @name Expression traits.
    /// @{

    /// Validates if the active type results of two expressions are the same or compatible. `void` results are
    /// interpreted as the active type result of a constant expression.
    template<typename ResultA, typename ResultB, typename = void>
    struct ValidateResultImpl {
      private:
        static bool constexpr isAVoid = std::is_same<void, ResultA>::value;
        static bool constexpr isBVoid = std::is_same<void, ResultB>::value;
        static bool constexpr isBothVoid = isAVoid & isBVoid;
        static bool constexpr isBothSame = std::is_same<ResultA, ResultB>::value;

        // Either one can be void (aka. constant value) but not both otherwise both need to be the same.
        CODI_STATIC_ASSERT((!isBothVoid) & (!isAVoid | !isBVoid | isBothSame), "Result types need to be the same.");

      public:

        /// The resulting active type of an expression.
        using ActiveResult = typename std::conditional<isBVoid, ResultA, ResultB>::type;
    };

    /// \copydoc ValidateResultImpl
    template<typename ResultA, typename ResultB>
    using ValidateResult = ValidateResultImpl<ResultA, ResultB>;

    /// @}
    /*******************************************************************************/
    /// @name Detection of specific node types
    /// @{

    /// If the expression inherits from ExpressionInterface. Is either std::false_type or std::true_type
    template<typename Expr, typename = void>
    struct IsExpression : std::false_type {};

#ifndef DOXYGEN_DISABLE
    template<typename Expr>
    struct IsExpression<Expr, typename enable_if_base_of<ExpressionInterface<typename Expr::Real, Expr>, Expr>::type>
        : std::true_type {};
#endif

#if CODI_IS_CPP14
    /// Value entry of IsExpression
    template<typename Expr>
    bool constexpr isExpression = IsExpression<Expr>::value;
#endif

    /// Enable if wrapper for IsExpression
    template<typename Expr, typename T = void>
    using EnableIfExpression = typename std::enable_if<IsExpression<Expr>::value, T>::type;

    /// If the expression inherits from LhsExpressionInterface. Is either std::false_type or std::true_type
    template<typename Expr, typename = void>
    struct IsLhsExpression : std::false_type {};

#ifndef DOXYGEN_DISABLE
    template<typename Expr>
    struct IsLhsExpression<
        Expr, typename enable_if_base_of<
                  LhsExpressionInterface<typename Expr::Real, typename Expr::Gradient, typename Expr::Tape, Expr>,
                  Expr>::type> : std::true_type {};

    template<typename Tape>
    struct IsLhsExpression<StaticContextActiveType<Tape>> : std::true_type {};
#endif

#if CODI_IS_CPP14
    /// Value entry of IsLhsExpression
    template<typename Expr>
    bool constexpr isLhsExpression = IsLhsExpression<Expr>::value;
#endif

    /// Enable if wrapper for IsLhsExpression
    template<typename Expr, typename T = void>
    using EnableIfLhsExpression = typename std::enable_if<IsLhsExpression<Expr>::value, T>::type;

    /// If the expression inherits from ConstantExpression. Is either std::false_type or std::true_type
    template<typename Expr>
    struct IsConstantExpression : std::false_type {};

#ifndef DOXYGEN_DISABLE
    template<typename Real, template<typename> class ConversionOperator>
    struct IsConstantExpression<ConstantExpression<Real, ConversionOperator>> : std::true_type {};
#endif

#if CODI_IS_CPP14
    template<typename Expr>
    /// Value entry of IsConstantExpression
    bool constexpr isConstantExpression = IsConstantExpression<Expr>::value;
#endif

    /// Enable if wrapper for IsConstantExpression
    template<typename Expr, typename T = void>
    using EnableIfConstantExpression = typename std::enable_if<IsConstantExpression<Expr>::value, T>::type;

    /// If the expression inherits from StaticContextActiveType. Is either std::false_type or std::true_type
    template<typename Expr>
    struct IsStaticContextActiveType : std::false_type {};

#ifndef DOXYGEN_DISABLE
    template<typename Tape>
    struct IsStaticContextActiveType<StaticContextActiveType<Tape>> : std::true_type {};
#endif

#if CODI_IS_CPP14
    /// Value entry of IsStaticContextActiveType
    template<typename Expr>
    bool constexpr isStaticContextActiveType = IsStaticContextActiveType<Expr>::value;
#endif

    /// Enable if wrapper for IsStaticContextActiveType
    template<typename Expr, typename T = void>
    using EnableIfStaticContextActiveType = typename std::enable_if<IsStaticContextActiveType<Expr>::value, T>::type;

    /// @}
    /*******************************************************************************/
    /// @name Static values on expressions
    /// @{

    /// Counts the number of nodes that inherit from LhsExpressionInterface in the expression.
    template<typename Expr>
    struct NumberOfActiveTypeArguments : public CompileTimeTraversalLogic<size_t, NumberOfActiveTypeArguments<Expr>> {
      public:

        /// Abbreviation for the base class type.
        using Base = CompileTimeTraversalLogic<size_t, NumberOfActiveTypeArguments<Expr>>;

        /// \copydoc CompileTimeTraversalLogic::leaf()
        template<typename Node, typename = ExpressionTraits::EnableIfLhsExpression<Node>>
        CODI_INLINE static size_t constexpr leaf() {
          return 1;
        }
        using CompileTimeTraversalLogic<size_t, NumberOfActiveTypeArguments>::leaf;

        /// See NumberOfActiveTypeArguments
        static size_t constexpr value = Base::template eval<Expr>();
    };

#if CODI_IS_CPP14
    /// Value entry of NumberOfActiveTypeArguments
    template<typename Expr>
    bool constexpr numberOfActiveTypeArguments = NumberOfActiveTypeArguments<Expr>::value;
#endif

    /// Counts the number of types that inherit from ConstantExpression in the expression.
    template<typename Expr>
    struct NumberOfConstantTypeArguments
        : public CompileTimeTraversalLogic<size_t, NumberOfConstantTypeArguments<Expr>> {
      public:

        /// Abbreviation for the base class type.
        using Base = CompileTimeTraversalLogic<size_t, NumberOfConstantTypeArguments<Expr>>;

        /// \copydoc CompileTimeTraversalLogic::leaf()
        template<typename Node, typename = EnableIfConstantExpression<Node>>
        CODI_INLINE static size_t constexpr leaf() {
          return 1;
        }
        using CompileTimeTraversalLogic<size_t, NumberOfConstantTypeArguments>::leaf;

        /// See NumberOfConstantTypeArguments
        static size_t constexpr value = Base::template eval<Expr>();
    };

#if CODI_IS_CPP14
    /// Value entry of NumberOfConstantTypeArguments
    template<typename Expr>
    bool constexpr numberOfConstantTypeArguments = NumberOfConstantTypeArguments<Expr>::value;
#endif

    /// Counts the number of nodes in the expression.
    template<typename Expr>
    struct NumberOfOperations : public CompileTimeTraversalLogic<size_t, NumberOfOperations<Expr>> {
      public:

        /// \copydoc CompileTimeTraversalLogic::node()
        template<typename Node>
        CODI_INLINE static size_t constexpr node() {
          return 1 + NumberOfOperations::template toLinks<Node>();
        }

        /// See NumberOfOperations
        static size_t constexpr value = NumberOfOperations::template eval<Expr>();
    };

#if CODI_IS_CPP14
    /// Value entry of NumberOfOperations
    template<typename Expr>
    bool constexpr numberOfOperations = NumberOfOperations<Expr>::value;
#endif

    /// @}
  }
}
