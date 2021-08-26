#pragma once

#include <iostream>

#include "../aux/macros.hpp"
#include "../config.h"
#include "../traits/expressionTraits.hpp"
#include "../traits/realTraits.hpp"
#include "logic/nodeInterface.hpp"

/** \copydoc codi::Namespace */
namespace codi {

  /**
   * @brief Base class for all CoDiPack expressions.
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

      using Real = CODI_DD(_Real, double);               ///< See ExpressionInterface.
      using Impl = CODI_DD(_Impl, ExpressionInterface);  ///< See ExpressionInterface.

      using ActiveResult = CODI_UNDEFINED;  ///< Type into which the expression can be converted. Usually also the type
                                            ///< from which is constructed.

      /// Cast to the implementation.
      CODI_INLINE Impl const& cast() const {
        return static_cast<Impl const&>(*this);
      }

#if CODI_ImplicitConversion
      /// Implicit cast for CoDiPack expressions.
      CODI_INLINE operator const Real() const {
        Warning::implicitCast<Config::ImplicitConversionWarning>();

        return cast().getValue();
      }
#endif

      /*******************************************************************************/
      /// @name Interface definition
      /// @{

      using StoreAs = ExpressionInterface;  ///< Defines how this expression is stored in an expression tree.

      /// Compute the primal value that is usually evaluated by the statement/expression.
      CODI_INLINE Real const getValue() const;

      /// Get the Jacobian with respect to the given argument.
      ///
      /// This is just the local Jacobian and not the one for the whole expression tree.
      template<size_t argNumber>
      CODI_INLINE Real getJacobian() const;

      /// @}

    private:
      ExpressionInterface& operator=(ExpressionInterface const&) = delete;
  };

#ifndef DOXYGEN_DISABLE
  template<typename _Type>
  struct RealTraits::TraitsImplementation<_Type, ExpressionTraits::EnableIfExpression<_Type>> {
    public:

      using Type = CODI_DD(_Type, CODI_T(ExpressionInterface<double, _Type>));
      using Real = typename Type::Real;

      using PassiveReal = RealTraits::PassiveReal<Real>;

      static int constexpr MaxDerivativeOrder = 1 + RealTraits::MaxDerivativeOrder<Real>();

      static CODI_INLINE PassiveReal const& getPassiveValue(Type const& v) {
        return RealTraits::getPassiveValue(v.getValue());
      }
  };
#endif

  /// Write the primal value to the stream.
  template<typename Expr>
  ExpressionTraits::EnableIfExpression<Expr, std::ostream>& operator<<(std::ostream& out, Expr const& v) {
    out << v.getValue();

    return out;
  }
}
