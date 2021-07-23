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
   * @tparam _Real Type of the Jacobian values.
   * @tparam _Impl  Class implementing this interface.
   */
  template<typename _Real, typename _Impl>
  struct JacobianComputationLogic : public TraversalLogic<_Impl> {
    public:

      using Real = CODI_DD(_Real, double);                            ///< See JacobianComputationLogic.
      using Impl = CODI_DD(_Impl, CODI_T(TraversalLogic<CODI_ANY>));  ///< See JacobianComputationLogic.

      /*******************************************************************************/
      /// @name Interface definition
      /// @{

      /// Called for leaf nodes which implement LhsExpressionInterface.
      template<typename Node, typename... Args>
      void handleJacobianOnActive(Node const& node, Real jacobian, Args&&... args);

      /// @}
      /*******************************************************************************/
      /// @name Overwrites from TraversalLogic
      /// @{

      /// \copydoc codi::TraversalLogic::leaf()
      template<typename Node, typename... Args>
      CODI_INLINE ExpressionTraits::EnableIfLhsExpression<Node> leaf(Node const& node, Real jacobian, Args&&... args) {
        cast().handleJacobianOnActive(node, jacobian, std::forward<Args>(args)...);
      }

      using TraversalLogic<Impl>::leaf;

      /// Computes the \ref sec_reverseAD "reverse" AD equation for this link.
      ///
      /// The Jacobian is multiplied with the Jacobian of the link. The result is forwarded to the child.
      template<size_t ChildNumber, typename Child, typename Root, typename... Args>
      CODI_INLINE void link(Child const& child, Root const& root, Real const& jacobian, Args&&... args) {
        Real curJacobian = root.template getJacobian<ChildNumber>() * jacobian;

        cast().toNode(child, curJacobian, std::forward<Args>(args)...);
      }

      /// @}

    private:

      CODI_INLINE Impl& cast() {
        return static_cast<Impl&>(*this);
      }
  };
}
