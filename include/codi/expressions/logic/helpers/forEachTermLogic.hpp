#pragma once

#include <type_traits>
#include <utility>

#include "../../../aux/macros.hpp"
#include "../../../config.h"
#include "../../../traits/expressionTraits.hpp"
#include "../traversalLogic.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Implement logic for termination nodes only.
   *
   * This class calls:
   *  - handleActive for each termination node that implements LhsExpressionInterface
   *  - handleConstant for each termination node that implements ConstantExpression
   *
   * For details about the expression traversal see TraversalLogic.
   *
   * @tparam _Impl  Class implementing this interface.
   */
  template<typename _Impl>
  struct ForEachTermLogic : public TraversalLogic<_Impl> {
    public:

      using Impl = CODI_DD(_Impl, CODI_T(TraversalLogic<CODI_ANY>));  ///< See ForEachTermLogic

      /*******************************************************************************/
      /// @name Interface definition
      /// @{

      /// Called for termination nodes which implement LhsExpressionInterface
      template<typename Node, typename... Args>
      void handleActive(Node const& node, Args&&... args) {
        CODI_UNUSED(node, args...);
      }

      /// Called for termination nodes which implement ConstantExpression
      template<typename Node, typename... Args>
      void handleConstant(Node const& node, Args&&... args) {
        CODI_UNUSED(node, args...);
      }

      /// @}
      /*******************************************************************************/
      /// @name Overwrites from TraversalLogic
      /// @{

      /// \copydoc codi::TraversalLogic::term()
      template<typename Node, typename... Args>
      CODI_INLINE ExpressionTraits::EnableIfLhsExpression<Node> term(Node const& node, Args&&... args) {
        cast().handleActive(node, std::forward<Args>(args)...);
      }

      /// \copydoc codi::TraversalLogic::term()
      template<typename Node, typename... Args>
      CODI_INLINE ExpressionTraits::EnableIfConstantExpression<Node> term(Node const& node, Args&&... args) {
        cast().handleConstant(node, std::forward<Args>(args)...);
      }

      using TraversalLogic<Impl>::term;

    /// @}

    private:

      CODI_INLINE Impl& cast() {
        return static_cast<Impl&>(*this);
      }
  };
}
