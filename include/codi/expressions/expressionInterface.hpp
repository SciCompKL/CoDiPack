#pragma once

#include "../aux/macros.hpp"
#include "../config.h"
#include "logic/nodeInterface.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Base class for all CoDiPack expression.
   *
   * See \ref Expressions "Expression" design documentation for details about the expression system in CoDiPack.
   *
   * This interface resembles an rvalue in C++.
   *
   * @tparam _Real  Original primal value of the statement/expression.
   * @tparam _Impl  Class implementing this interface.
   */
  template<typename _Real, typename _Impl>
  struct ExpressionInterface : public NodeInterface<_Impl> {
    public:

      using Real = CODI_DD(_Real, double);  ///< See ExpressionInterface
      using Impl = CODI_DD(_Impl, ExpressionInterface);  ///< See ExpressionInterface

      /// Cast to the implementation.
      CODI_INLINE Impl const& cast() const {
        return static_cast<Impl const&>(*this);
      }

      /*******************************************************************************/
      /// @name Interface definition
      /// @{

      using StoreAs = ExpressionInterface;  ///< Defines how this expression is stored in an expression tree.

      /// Compute the primal value that is usually evaluated by the statement/expression.
      CODI_INLINE Real const getValue() const;

      /// Get the Jacobian with respect to the given argument.
      ///
      /// This is just the local Jacobian and not the one from the whole expression tree.
      template<size_t argNumber>
      CODI_INLINE Real getJacobian() const;

      /// @}

    private:
      ExpressionInterface& operator=(ExpressionInterface const&) = delete;
  };
}
