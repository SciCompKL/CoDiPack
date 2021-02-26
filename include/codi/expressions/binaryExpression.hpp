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
   * @brief Implements the logic for a BinaryExpression.
   *
   * See \ref Expressions "Expression" design documentation for details about the expression system in CoDiPack.
   *
   * @tparam _Real  Original primal value of the statement/expression.
   */
  template<typename _Real>
  struct BinaryOperation {
    public:

      using Real = CODI_DD(_Real, double);  ///< See BinaryOperation

      /// Compute the primal value from the arguments.
      ///
      /// The type of the arguments is the result of a getValue call on the expressions.
      template<typename ArgA, typename ArgB>
      static CODI_INLINE Real primal(ArgA const& argA, ArgB const& argB);

      /// Compute the gradient with respect to the first argument
      ///
      /// The type of the arguments is the result of a getValue call on the expressions.
      template<typename ArgA, typename ArgB>
      static CODI_INLINE Real gradientA(ArgA const& argA, ArgB const& argB, Real const& result);

      /// Compute the gradient with respect to the second argument
      ///
      /// The type of the arguments is the result of a getValue call on the expressions.
      template<typename ArgA, typename ArgB>
      static CODI_INLINE Real gradientB(ArgA const& argA, ArgB const& argB, Real const& result);
  };


  /**
   * @brief Represents an operator with two arguments in the expression tree.
   *
   * See \ref Expressions "Expression" design documentation for details about the expression system in CoDiPack.
   *
   * @tparam _Real  Original primal value of the statement/expression.
   * @tparam _ArgA  The ExpressionInterface type of the first argument.
   * @tparam _ArgB  The ExpressionInterface type of the second argument.
   * @tparam _Operation  The logic for computing the primal value and Jacobians. Needs to implement BinaryOperation.
   */
  template<typename _Real, typename _ArgA, typename _ArgB, template<typename> class _Operation>
  struct BinaryExpression : public ExpressionInterface<_Real, BinaryExpression<_Real, _ArgA, _ArgB, _Operation> > {
    public:
      using Real = CODI_DD(_Real, double);  ///< See BinaryExpression
      using ArgA = CODI_DD(_ArgA, CODI_T(ExpressionInterface<double, CODI_ANY>));  ///< See BinaryExpression
      using ArgB = CODI_DD(_ArgB, CODI_T(ExpressionInterface<double, CODI_ANY>));  ///< See BinaryExpression
      using Operation = CODI_DD(CODI_T(_Operation<Real>), CODI_T(BinaryOperation<Real>));  ///< See BinaryExpression

      typename ArgA::StoreAs argA;  ///< First argument of the expression
      typename ArgB::StoreAs argB;  ///< Second argument of the expression
      Real result;  ///< Precomputed result

      /// Constructor
      template<typename RealA, typename RealB>
      explicit BinaryExpression(ExpressionInterface<RealA, ArgA> const& argA, ExpressionInterface<RealB, ArgB> const& argB) :
        argA(argA.cast()),
        argB(argB.cast()),
        result(Operation::primal(this->argA.getValue(), this->argB.getValue())) {}

      /*******************************************************************************/
      /// @name Implementation of ExpressionInterface
      /// @{

      using StoreAs = BinaryExpression;  ///< \copydoc codi::ExpressionInterface::StoreAs

      /// \copydoc codi::ExpressionInterface::getValue()
      CODI_INLINE Real const& getValue() const {
        return result;
      }

      /// \copydoc codi::ExpressionInterface::getJacobian()
      template<size_t argNumber>
      CODI_INLINE Real getJacobian() const {
        if (0 == argNumber) {
          return Operation::gradientA(argA.getValue(), argB.getValue(), result);
        } else {
          return Operation::gradientB(argA.getValue(), argB.getValue(), result);
        }
      }

      /// @}
      /*******************************************************************************/
      /// @name Implementation of NodeInterface
      /// @{

      static bool constexpr EndPoint = false;  ///< \copydoc codi::NodeInterface::EndPoint

      /// \copydoc codi::NodeInterface::forEachLink
      template<typename Logic, typename... Args>
      CODI_INLINE void forEachLink(TraversalLogic<Logic>& logic, Args&&... args) const {
        logic.cast().template link<0>(argA, *this, std::forward<Args>(args)...);
        logic.cast().template link<1>(argB, *this, std::forward<Args>(args)...);
      }

      /// \copydoc codi::NodeInterface::forEachLinkConstExpr
      template<typename CompileTimeLogic, typename... Args>
      CODI_INLINE static typename CompileTimeLogic::ResultType constexpr forEachLinkConstExpr(Args&&... args) {
        return CompileTimeLogic::reduce(
              CompileTimeLogic::template link<0, ArgA, BinaryExpression>(std::forward<Args>(args)...),
              CompileTimeLogic::template link<1, ArgB, BinaryExpression>(std::forward<Args>(args)...));
      }

      /// @}
  };
}
