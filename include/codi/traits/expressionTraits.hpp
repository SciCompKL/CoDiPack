#pragma once

#include <type_traits>

#include "../aux/macros.h"
#include "../config.h"
#include "../expressions/logic/compileTimeTraversalLogic.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  template<typename _Real, typename _Gradient, typename _Tape, typename _Impl>
  struct LhsExpressionInterface;

  template<typename _Real>
  struct ConstantExpression;

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
    typename std::enable_if<
      std::is_base_of<
        LhsExpressionInterface<typename Impl::Real, typename Impl::Gradient, typename Impl::Tape, Impl>,
        Impl
      >::value
    >::type
  > : std::true_type {};


  template<typename Impl>
  using isLhsExpression = IsLhsExpression<Impl>;

  template<typename Impl>
  using enableIfLhsExpression = typename std::enable_if<isLhsExpression<Impl>::value>::type;

  template<typename Impl>
  struct IsConstantExpression : std::false_type {};

  template<typename Real>
  struct IsConstantExpression<ConstantExpression<Real>> : std::true_type {};

  template<typename Impl>
  using isConstantExpression = IsConstantExpression<Impl>;

  template<typename Impl>
  using enableIfConstantExpression = typename std::enable_if<isConstantExpression<Impl>::value>::type;

  /*******************************************************************************
   * Section: Static values on expressions
   *
   * Description: TODO
   *
   */

  template<typename Expr>
  struct MaxNumberOfActiveArguments : public CompileTimeTraversalLogic<size_t, MaxNumberOfActiveArguments<Expr>> {
    public:

      static size_t constexpr value = MaxNumberOfActiveArguments::template eval<Expr>();

      template<typename Node, typename = enableIfLhsExpression<Node>>
      CODI_INLINE static constexpr size_t term() {
        return 1;
      }
      using CompileTimeTraversalLogic<size_t, MaxNumberOfActiveArguments>::term;
  };

  template<typename Expr>
  struct MaxNumberOfConstantArguments : public CompileTimeTraversalLogic<size_t, MaxNumberOfConstantArguments<Expr>> {
    public:

      static size_t constexpr value = MaxNumberOfConstantArguments::template eval<Expr>();

      template<typename Node, typename = enableIfConstantExpression<Node>>
      CODI_INLINE static constexpr size_t term() {
        return 1;
      }
      using CompileTimeTraversalLogic<size_t, MaxNumberOfConstantArguments>::term;
  };
}
