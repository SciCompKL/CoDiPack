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

    /*******************************************************************************
     * Section: Detection of specific node types
     *
     * Description: TODO
     *
     */

    template<typename Impl, typename = void>
    struct IsLhsExpression : std::false_type {};

    template<typename Impl>
    struct IsLhsExpression<
      Impl,
      typename enable_if_base_of<
        LhsExpressionInterface<typename Impl::Real, typename Impl::Gradient, typename Impl::Tape, Impl>,
        Impl
      >::type
    > : std::true_type {};

    template<typename Tape>
    struct IsLhsExpression<StaticContextActiveType<Tape>> : std::true_type {};

#if CODI_IS_CPP14
    template<typename Impl>
    bool constexpr isLhsExpression = IsLhsExpression<Impl>::value;
#endif

    template<typename Impl>
    using EnableIfLhsExpression = typename std::enable_if<IsLhsExpression<Impl>::value>::type;

    template<typename Impl>
    struct IsConstantExpression : std::false_type {};

    template<typename Real>
    struct IsConstantExpression<ConstantExpression<Real>> : std::true_type {};

    template<typename Real, size_t offset>
    struct IsConstantExpression<StaticContextConstantExpression<Real, offset>> : std::true_type {};

#if CODI_IS_CPP14
    template<typename Impl>
    bool constexpr isConstantExpression = IsConstantExpression<Impl>::value;
#endif

    template<typename Impl>
    using EnableIfConstantExpression = typename std::enable_if<IsConstantExpression<Impl>::value>::type;

    template<typename Impl>
    struct IsStaticContextActiveType : std::false_type {};

    template<typename Tape>
    struct IsStaticContextActiveType<StaticContextActiveType<Tape>> : std::true_type {};

#if CODI_IS_CPP14
    template<typename Impl>
    bool constexpr isStaticContextActiveType = IsStaticContextActiveType<Impl>::value;
#endif

    template<typename Impl>
    using EnableIfStaticContextActiveType = typename std::enable_if<IsStaticContextActiveType<Impl>::value>::type;

    /*******************************************************************************
     * Section: Static values on expressions
     *
     * Description: TODO
     *
     */

    template<typename Expr>
    struct NumberOfActiveTypeArguments : public CompileTimeTraversalLogic<size_t, NumberOfActiveTypeArguments<Expr>> {
      public:

        template<typename Node, typename = ExpressionTraits::EnableIfLhsExpression<Node>>
        CODI_INLINE static size_t constexpr term() {
          return 1;
        }
        using CompileTimeTraversalLogic<size_t, NumberOfActiveTypeArguments>::term;

        static size_t constexpr value = NumberOfActiveTypeArguments::template eval<Expr>();
    };

#if CODI_IS_CPP14
    template<typename Expr>
    bool constexpr numberOfActiveTypeArguments = NumberOfActiveTypeArguments<Expr>::value;
#endif

    template<typename Expr>
    struct NumberOfConstantTypeArguments : public CompileTimeTraversalLogic<size_t, NumberOfConstantTypeArguments<Expr>> {
      public:

        template<typename Node, typename = EnableIfConstantExpression<Node>>
        CODI_INLINE static size_t constexpr term() {
          return 1;
        }
        using CompileTimeTraversalLogic<size_t, NumberOfConstantTypeArguments>::term;

        static size_t constexpr value = NumberOfConstantTypeArguments::template eval<Expr>();
    };

#if CODI_IS_CPP14
    template<typename Expr>
    bool constexpr numberOfConstantTypeArguments = NumberOfConstantTypeArguments<Expr>::value;
#endif

  }










}
