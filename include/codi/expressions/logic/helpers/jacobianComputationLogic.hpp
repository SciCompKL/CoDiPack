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
   *
   * @tparam T_Impl  Class implementing this interface.
   */
  template<typename T_Impl>
  struct JacobianComputationLogic : public TraversalLogic<T_Impl> {
    public:

      using Impl = CODI_DD(T_Impl, CODI_T(TraversalLogic<CODI_ANY>));  ///< See JacobianComputationLogic.

      /*******************************************************************************/
      /// @name Interface definition
      /// @{

      /// Called for leaf nodes which implement LhsExpressionInterface.
      template<typename Node, typename Jacobian, typename... Args>
      void handleJacobianOnActive(Node const& node, Jacobian jacobian, Args&&... args);

      /// @}
      /*******************************************************************************/
      /// @name Overwrites from TraversalLogic
      /// @{

      /// \copydoc codi::TraversalLogic::leaf()
      template<typename Node, typename Jacobian, typename... Args>
      CODI_INLINE ExpressionTraits::EnableIfLhsExpression<Node> leaf(Node const& node, Jacobian jacobian,
                                                                     Args&&... args) {
        cast().handleJacobianOnActive(node, jacobian, std::forward<Args>(args)...);
      }

      using TraversalLogic<Impl>::leaf;

      /// Computes the \ref sec_reverseAD "reverse" AD equation for this link.
      ///
      /// The Jacobian is multiplied with the Jacobian of the link. The result is forwarded to the child.
      template<size_t ChildNumber, typename Child, typename Root, typename Jacobian, typename... Args>
      CODI_INLINE void link(Child const& child, Root const& root, Jacobian const& jacobian, Args&&... args) {
        auto linkJacobian = root.template getJacobian<ChildNumber>();
        auto curJacobian = linkJacobian * jacobian;

        cast().toNode(child, curJacobian, std::forward<Args>(args)...);
      }

      /// @}

    private:

      CODI_INLINE Impl& cast() {
        return static_cast<Impl&>(*this);
      }
  };
}
