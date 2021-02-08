#pragma once

#include "../aux/macros.hpp"
#include "../config.h"
#include "expressionInterface.hpp"
#include "logic/compileTimeTraversalLogic.hpp"
#include "logic/nodeInterface.hpp"
#include "logic/traversalLogic.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Implements the logic for a UnaryExpression.
   *
   * See \ref Expressions "Expression" design documentation for details about the expression system in CoDiPack.
   *
   * @tparam _Real  Original primal value of the statement/expression.
   */
  template<typename _Real>
  struct UnaryOperation {
    public:

      using Real = CODI_DECLARE_DEFAULT(_Real, double); ///< See UnaryOperation

      /// Compute the primal value from the argument.
      ///
      /// The type of the argument is the result of a getValue call on the expression.
      template<typename Arg>
      static CODI_INLINE Real primal(Arg const& arg);

      /// Compute the gradient with respect to the argument
      ///
      /// The type of the argument is the result of a getValue call on the expression.
      template<typename Arg>
      static CODI_INLINE Real gradient(Arg const& arg, Real const& result);
  };

  /**
   * @brief Represents an operator with one argument in the expression tree.
   *
   * See \ref Expressions "Expression" design documentation for details about the expression system in CoDiPack.
   *
   * @tparam _Real  Original primal value of the statement/expression.
   * @tparam _Arg  The ExpressionInterface type of the argument.
   * @tparam _Operation  The logic for computing the primal value and Jacobian. Needs to implement UnaryOperation.
   */
  template<typename _Real, typename _Arg, template<typename> class _Operation>
  struct UnaryExpression : public ExpressionInterface<_Real, UnaryExpression<_Real, _Arg, _Operation> > {
    public:

      using Real = CODI_DECLARE_DEFAULT(_Real, double); ///< See UnaryExpression
      using Arg = CODI_DECLARE_DEFAULT(_Arg, CODI_TEMPLATE(ExpressionInterface<double, CODI_ANY>)); ///< See UnaryExpression
      using Operation = CODI_DECLARE_DEFAULT(CODI_TEMPLATE(_Operation<Real>), CODI_TEMPLATE(UnaryOperation<Real>)); ///< See UnaryExpression

      typename Arg::StoreAs arg; ///< Argument of the expression
      Real result; ///< Precomputed result

      /// Constructor
      template<typename RealArg>
      explicit UnaryExpression(ExpressionInterface<RealArg, Arg> const& arg) :
        arg(arg.cast()),
        result(Operation::primal(this->arg.getValue())) {}

      /*******************************************************************************/
      /// @name Implementation of ExpressionInterface
      /// @{

      using StoreAs = UnaryExpression; ///< \copydoc codi::ExpressionInterface::StoreAs

      /// \copydoc codi::ExpressionInterface::getValue()
      CODI_INLINE Real const& getValue() const {
        return result;
      }

      /// \copydoc codi::ExpressionInterface::getJacobian()
      template<size_t argNumber>
      CODI_INLINE Real getJacobian() const {
        return Operation::gradient(arg.getValue(), result);
      }

      /// @}
      /*******************************************************************************/
      /// @name Implementation of NodeInterface
      /// @{

      static bool constexpr EndPoint = false; ///< \copydoc codi::NodeInterface::EndPoint

      /// \copydoc codi::NodeInterface::forEachLink
      template<typename Logic, typename ... Args>
      CODI_INLINE void forEachLink(TraversalLogic<Logic>& logic, Args&& ... args) const {
        logic.cast().template link<0>(arg, *this, std::forward<Args>(args)...);
      }

      /// \copydoc codi::NodeInterface::forEachLinkConstExpr
      template<typename CompileTimeLogic, typename ... Args>
      CODI_INLINE static typename CompileTimeLogic::ResultType constexpr forEachLinkConstExpr(Args&& ... args) {
        return CompileTimeLogic::template link<0, Arg, UnaryExpression>(std::forward<Args>(args)...);
      }

      /// @}
  };
}
