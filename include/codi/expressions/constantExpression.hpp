#pragma once

#include "../aux/macros.hpp"
#include "../config.h"
#include "expressionInterface.hpp"
#include "logic/nodeInterface.hpp"
#include "logic/traversalLogic.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Represents constant values in the expression tree.
   *
   * All values that are not a CoDiPack type are considered constant, for example values like 4.0 or double a.
   *
   * See \ref Expressions "Expression" design documentation for details about the expression system in CoDiPack.
   *
   * @tparam T_Real  Original primal value of the statement/expression.
   */
  template<typename T_Real>
  struct ConstantExpression : public ExpressionInterface<T_Real, ConstantExpression<T_Real>> {
    public:

      using Real = CODI_DD(T_Real, double);  ///< See ConstantExpression.

    private:
      Real primalValue;

    public:

      /// Constructor
      CODI_INLINE ConstantExpression(Real const& v) : primalValue(v) {}

      /*******************************************************************************/
      /// @name Implementation of ExpressionInterface
      /// @{

      using StoreAs = ConstantExpression;  ///< \copydoc codi::ExpressionInterface::StoreAs
      using ActiveResult = void;           ///< \copydoc codi::ExpressionInterface::ActiveResult

      /// \copydoc codi::ExpressionInterface::getValue
      CODI_INLINE Real const& getValue() const {
        return primalValue;
      }

      /// \copydoc codi::ExpressionInterface::getJacobian()
      template<size_t argNumber>
      CODI_INLINE Real getJacobian() const {
        return Real();
      }

      /// @}
      /*******************************************************************************/
      /// @name Implementation of NodeInterface
      /// @{

      static bool constexpr EndPoint = true;  ///< \copydoc codi::NodeInterface::EndPoint

      /// \copydoc codi::NodeInterface::forEachLink()
      template<typename Logic, typename... Args>
      CODI_INLINE void forEachLink(TraversalLogic<Logic>& logic, Args&&... args) const {
        CODI_UNUSED(logic, args...);
      }

      /// \copydoc codi::NodeInterface::forEachLinkConstExpr()
      template<typename Logic, typename... Args>
      CODI_INLINE static typename Logic::ResultType constexpr forEachLinkConstExpr(Args&&... CODI_UNUSED_ARG(args)) {
        return Logic::NeutralElement;
      }

      /// @}
  };
}
