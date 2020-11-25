#pragma once

#include <type_traits>

#include "../aux/macros.hpp"
#include "../config.h"
#include "../expressions/logic/compileTimeTraversalLogic.hpp"
#include "aux/enableIfHelpers.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  template<typename _Real, typename _Gradient, typename _Tape, typename _Impl>
  struct LhsExpressionInterface;

  template<typename _Real>
  struct ConstantExpression;

  template<typename _Tape>
  struct StaticContextActiveType;

  template<typename _Real, size_t _offset>
  struct StaticContextConstantExpression;

  namespace ExpressionTraits {

    /*******************************************************************************/
    /// @name Detection of specific node types
    /// @{

    /// If the expression inherits from LhsExpressionInterface. Is either std::false_type or std::true_type
    template<typename Expr, typename = void>
    struct IsLhsExpression : std::false_type {};

#ifndef DOXYGEN_DISABLE
    template<typename Expr>
    struct IsLhsExpression<
      Expr,
      typename enable_if_base_of<
        LhsExpressionInterface<typename Expr::Real, typename Expr::Gradient, typename Expr::Tape, Expr>,
        Expr
      >::type
    > : std::true_type {};

    template<typename Tape>
    struct IsLhsExpression<StaticContextActiveType<Tape>> : std::true_type {};
#endif

#if CODI_IS_CPP14
    /// Value entry of IsLhsExpression
    template<typename Expr>
    bool constexpr isLhsExpression = IsLhsExpression<Expr>::value;
#endif

    /// Enable if wrapper for IsLhsExpression
    template<typename Expr>
    using EnableIfLhsExpression = typename std::enable_if<IsLhsExpression<Expr>::value>::type;

    /// If the expression inherits from ConstantExpression. Is either std::false_type or std::true_type
    template<typename Expr>
    struct IsConstantExpression : std::false_type {};

#ifndef DOXYGEN_DISABLE
    template<typename Real>
    struct IsConstantExpression<ConstantExpression<Real>> : std::true_type {};

    template<typename Real, size_t offset>
    struct IsConstantExpression<StaticContextConstantExpression<Real, offset>> : std::true_type {};
#endif

#if CODI_IS_CPP14
    template<typename Expr>
    /// Value entry of IsConstantExpression
    bool constexpr isConstantExpression = IsConstantExpression<Expr>::value;
#endif

    /// Enable if wrapper for IsConstantExpression
    template<typename Expr>
    using EnableIfConstantExpression = typename std::enable_if<IsConstantExpression<Expr>::value>::type;

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
    template<typename Expr>
    using EnableIfStaticContextActiveType = typename std::enable_if<IsStaticContextActiveType<Expr>::value>::type;

    /// @}
    /*******************************************************************************/
    /// @name Static values on expressions
    /// @{

    /// Counts the number of nodes that inherit from LhsExpressionInterface in the expression.
    template<typename Expr>
    struct NumberOfActiveTypeArguments : public CompileTimeTraversalLogic<size_t, NumberOfActiveTypeArguments<Expr>> {
      public:

        /// \copydoc CompileTimeTraversalLogic::term()
        template<typename Node, typename = ExpressionTraits::EnableIfLhsExpression<Node>>
        CODI_INLINE static size_t constexpr term() {
          return 1;
        }
        using CompileTimeTraversalLogic<size_t, NumberOfActiveTypeArguments>::term;

        /// See NumberOfActiveTypeArguments
        static size_t constexpr value = NumberOfActiveTypeArguments::template eval<Expr>();
    };

#if CODI_IS_CPP14
    /// Value entry of NumberOfActiveTypeArguments
    template<typename Expr>
    bool constexpr numberOfActiveTypeArguments = NumberOfActiveTypeArguments<Expr>::value;
#endif

    /// Counts the number of types that inherit from ConstantExpression in the expression.
    template<typename Expr>
    struct NumberOfConstantTypeArguments : public CompileTimeTraversalLogic<size_t, NumberOfConstantTypeArguments<Expr>> {
      public:

        /// \copydoc CompileTimeTraversalLogic::term()
        template<typename Node, typename = EnableIfConstantExpression<Node>>
        CODI_INLINE static size_t constexpr term() {
          return 1;
        }
        using CompileTimeTraversalLogic<size_t, NumberOfConstantTypeArguments>::term;

        /// See NumberOfConstantTypeArguments
        static size_t constexpr value = NumberOfConstantTypeArguments::template eval<Expr>();
    };

#if CODI_IS_CPP14
    /// Value entry of NumberOfConstantTypeArguments
    template<typename Expr>
    bool constexpr numberOfConstantTypeArguments = NumberOfConstantTypeArguments<Expr>::value;
#endif

    /// @}
  }
}
