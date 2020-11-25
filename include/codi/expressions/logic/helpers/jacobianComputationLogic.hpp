#pragma once

#include <type_traits>
#include <utility>

#include "../../../aux/macros.hpp"
#include "../../../config.h"
#include "../traversalLogic.hpp"
#include "../../../traits/expressionTraits.hpp"

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

      using Real = CODI_DECLARE_DEFAULT(_Real, double); ///< See JacobianComputationLogic
      using Impl = CODI_DECLARE_DEFAULT(_Impl, CODI_TEMPLATE(TraversalLogic<CODI_ANY>)); ///< See JacobianComputationLogic

      /*******************************************************************************/
      /// @name Interface definition
      /// @{

      /// Called for termination nodes which implement LhsExpressionInterface.
      template<typename Node, typename ... Args>
      void handleJacobianOnActive(Node const& node, Real jacobian, Args&& ... args);

      /// @}
      /*******************************************************************************/
      /// @name Overwrites from TraversalLogic
      /// @{

      /// \copydoc codi::TraversalLogic::term()
      template<typename Node, typename ... Args>
      CODI_INLINE ExpressionTraits::enableIfLhsExpression<Node> term(Node const& node, Real jacobian, Args&& ... args) {
        cast().handleJacobianOnActive(node, jacobian, std::forward<Args>(args)...);
      }

      using TraversalLogic<Impl>::term;

      /// Computes the \ref sec_reverseAD "reverse" AD equation for this link.
      ///
      /// The Jacobian is multiplied with the Jacobian of the link. The result is forwarded to the leaf.
      template<size_t LeafNumber, typename Leaf, typename Root, typename ... Args>
      CODI_INLINE void link(Leaf const& leaf, Root const& root, Real const& jacobian, Args&& ... args) {

        Real curJacobian = root.template getJacobian<LeafNumber>() * jacobian;

        cast().toNode(leaf, curJacobian, std::forward<Args>(args)...);
      }

      /// @}

    private:

      CODI_INLINE Impl& cast() {
        return static_cast<Impl&>(*this);
      }

  };
}
