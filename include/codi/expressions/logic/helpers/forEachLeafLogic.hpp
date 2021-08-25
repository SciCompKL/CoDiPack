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
   * @brief Implement logic for leaf nodes only.
   *
   * This class calls:
   *  - handleActive for each leaf node that implements LhsExpressionInterface,
   *  - handleConstant for each leaf node that implements ConstantExpression.
   *
   * For details about the expression traversal see TraversalLogic.
   *
   * @tparam _Impl  Class implementing this interface.
   */
  template<typename _Impl>
  struct ForEachLeafLogic : public TraversalLogic<_Impl> {
    public:

      using Impl = CODI_DD(_Impl, CODI_T(TraversalLogic<CODI_ANY>));  ///< See ForEachLeafLogic.

      /*******************************************************************************/
      /// @name Interface definition
      /// @{

      /// Called for leaf nodes which implement LhsExpressionInterface.
      template<typename Node, typename... Args>
      void handleActive(Node const& node, Args&&... args) {
        CODI_UNUSED(node, args...);
      }

      /// Called for leaf nodes which implement ConstantExpression.
      template<typename Node, typename... Args>
      void handleConstant(Node const& node, Args&&... args) {
        CODI_UNUSED(node, args...);
      }

      /// @}
      /*******************************************************************************/
      /// @name Overwrites from TraversalLogic
      /// @{

      /// \copydoc codi::TraversalLogic::leaf()
      template<typename Node, typename... Args>
      CODI_INLINE ExpressionTraits::EnableIfLhsExpression<Node> leaf(Node const& node, Args&&... args) {
        cast().handleActive(node, std::forward<Args>(args)...);
      }

      /// \copydoc codi::TraversalLogic::leaf()
      template<typename Node, typename... Args>
      CODI_INLINE ExpressionTraits::EnableIfConstantExpression<Node> leaf(Node const& node, Args&&... args) {
        cast().handleConstant(node, std::forward<Args>(args)...);
      }

      /// @}

    private:

      CODI_INLINE Impl& cast() {
        return static_cast<Impl&>(*this);
      }
  };
}
